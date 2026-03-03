module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <future>
#include <iostream>
#include <numbers>
#include <optional>
#include <span>
#include <thread>
#include <vector>
#include <fstream>
#include "bblock/system.h"
#endif

module interactions_mbx;

#ifndef USE_LEGACY_HEADERS
import <numbers>;
import <iostream>;
import <algorithm>;
import <vector>;
import <span>;
import <cmath>;
import <optional>;
import <thread>;
import <future>;
import <fstream>;
#include "bblock/system.h"
#endif

import energy_status;
import potential_energy_vdw;
import potential_gradient_vdw;
import potential_energy_coulomb;
import potential_gradient_coulomb;
import potential_correction_vdw;
import potential_electrostatics;
import simulationbox;
import double3;
import double3x3;
import atom;
import energy_factor;
import gradient_factor;
import energy_status_inter;
import running_energy;
import component;
import units;
import threadpool;
import system;
import component;
import framework;

// Used to calculate the total energy of the system for the given fixed positions of framework and adsorbate molecules
// TO-DO: consider skip this process if the system is just empty box
double Interactions::computeFrameworkElecPermMBXEnergy(
        const System &system,
        const SimulationBox &simulationBox,
        const std::optional<Framework> &framework,
        std::span<const Atom> frameworkAtoms
) noexcept
{
    /***
    * Setup MBX System object.
    ***/
    bblock::System* mbx = new bblock::System();
    int cumulativeTagIndex = 1;
    
    // no molecules are added since only calculate the empty framework E_elec
    int numAtoms{0};
    /***
    * MBX Miscellaneous settings.
    ***/
    // We now pass MBX file from RASPA.
    std::string json_path = system.mbxSettingsFilePath;
    std::ifstream t(json_path);
    t.seekg(0, std::ios::end);
    int size = t.tellg();
    std::string json_settings;
    json_settings.resize(size);
    t.seekg(0);
    t.read(&json_settings[0], size);
    mbx->SetUpFromJson(json_settings);
    /***
    * Adding framework atoms.
    ***/
    // Framework  --> In RASPA3 sometimes the framework is std::nullopt
    // Need to know MBX behaviour without the framework
    if (framework.has_value())
    {
        std::vector<double> frameworkCoords(frameworkAtoms.size() * 3);
        std::vector<double> frameworkCharges(frameworkAtoms.size());
        std::vector<size_t> frameworkIsLocals(frameworkAtoms.size(), 1);
        std::vector<int> frameworkTags(frameworkAtoms.size());

        for (int i = 0; i < frameworkAtoms.size(); i++) {
            const Atom &atom = frameworkAtoms[i];

            frameworkCoords[i * 3] = atom.position.x;
            frameworkCoords[i * 3 + 1] = atom.position.y;
            frameworkCoords[i * 3 + 2] = atom.position.z;

            frameworkCharges[i] = atom.charge;              
            frameworkTags[i] = cumulativeTagIndex;

            cumulativeTagIndex += 1;
        }
    
        mbx->SetExternalChargesAndPositions(frameworkCharges, frameworkCoords, frameworkIsLocals, frameworkTags);
    }

    /***
    * MBX simulation box settings.
    ***/
    std::vector<double> box(9, 0.0);

    // Box x vector, used to be m11, m12, m13
    box[0] = simulationBox.cell.mm[0][0];
    box[1] = simulationBox.cell.mm[0][1];
    box[2] = simulationBox.cell.mm[0][2];

    // Box y vector, used to be m21, m22, m23
    box[3] = simulationBox.cell.mm[1][0];
    box[4] = simulationBox.cell.mm[1][1];
    box[5] = simulationBox.cell.mm[1][2];

    // Box z vector, used to be m31, m32, m33
    box[6] = simulationBox.cell.mm[2][0];
    box[7] = simulationBox.cell.mm[2][1];
    box[8] = simulationBox.cell.mm[2][2];
    
    for (int i = 0; i < 9; i++) {
        if (std::abs(box[i]) < 1e-10) {
            box[i] = 0.0; // Avoid numerical issues with very small values
        }
    }

    mbx->SetPBC(box);

    /***
     * Calculate MBX energy.
    ***/

    bool do_grads = false;
    mbx->Electrostatics(do_grads);
    double eelec_perm_framework = mbx->GetPermanentElectrostaticEnergy();
    /***
     * Destroy MBX System object.
    ***/

    delete mbx;

    return eelec_perm_framework;
}


// Used to calculate the total energy of the system for the given fixed positions of framework and adsorbate molecules
RunningEnergy Interactions::computeMBXEnergySystem(
        const System &system,
        const std::vector<Component> &components,
        const SimulationBox &simulationBox,
        const std::optional<Framework> &framework,
        std::span<const Atom> frameworkAtoms,
        std::span<const Atom> moleculeAtoms
) noexcept
{
    /***
    * Setup MBX System object.
    ***/
    bblock::System* mbx = new bblock::System();
    int cumulativeTagIndex = 1;
    
    /***
    * Adding molecule monomers.
    ***/
    // Adding monomers molecules for all the components
    int numAtoms{0};
    for (size_t i = 0; i < components.size(); ++i)
    {
        numAtoms = static_cast<int>(components[i].atoms.size());           // Not sure whether MBX can use size_t
        std::vector<double> atomCoordinates(3 * numAtoms, 0.0);
        std::vector<std::string> atomNames(numAtoms, "none");
        std::string compName = components[i].name;

        // Index of first molecule of the given component.
        size_t idx_start = system.indexOfFirstMolecule(i);
        idx_start -= frameworkAtoms.size();                                 // Remove number of framework atoms
        
        // Looping over all the molecules of given components
        for (size_t j = 0; j < system.numberOfMoleculesPerComponent[i]; ++j)
        {    
            for(size_t k = 0; k < numAtoms; ++k) 
            {
                atomCoordinates[k * 3] = moleculeAtoms[idx_start + k + j*numAtoms].position.x;
                atomCoordinates[k * 3 + 1] = moleculeAtoms[idx_start + k + j*numAtoms].position.y;
                atomCoordinates[k * 3 + 2] = moleculeAtoms[idx_start + k + j*numAtoms].position.z;
                atomNames[k] = std::to_string(moleculeAtoms[idx_start + k + j*numAtoms].type);
            }

            size_t islocal = 1;

            size_t tagIndex = cumulativeTagIndex;
            
            mbx->AddMonomer(atomCoordinates, atomNames, compName, islocal, tagIndex);

            cumulativeTagIndex += numAtoms;
        }
    }
    /***
    * MBX Miscellaneous settings.
    ***/
    // We now pass MBX file from RASPA.
    std::string json_path = system.mbxSettingsFilePath;
    std::ifstream t(json_path);
    t.seekg(0, std::ios::end);
    int size = t.tellg();
    std::string json_settings;
    json_settings.resize(size);
    t.seekg(0);
    t.read(&json_settings[0], size);
    mbx->SetUpFromJson(json_settings);
    /***
    * Adding framework atoms.
    ***/
    // Framework  --> In RASPA3 sometimes the framework is std::nullopt
    // Need to know MBX behaviour without the framework
    double eelec_perm_framework = 0.0;
    if (framework.has_value())
    {
        std::vector<double> frameworkCoords(frameworkAtoms.size() * 3);
        std::vector<double> frameworkCharges(frameworkAtoms.size());
        std::vector<size_t> frameworkIsLocals(frameworkAtoms.size(), 1);
        std::vector<int> frameworkTags(frameworkAtoms.size());

        for (int i = 0; i < frameworkAtoms.size(); i++) {
            const Atom &atom = frameworkAtoms[i];

            frameworkCoords[i * 3] = atom.position.x;
            frameworkCoords[i * 3 + 1] = atom.position.y;
            frameworkCoords[i * 3 + 2] = atom.position.z;

            frameworkCharges[i] = atom.charge;              
            frameworkTags[i] = cumulativeTagIndex;

            cumulativeTagIndex += 1;
        }
    
        mbx->SetExternalChargesAndPositions(frameworkCharges, frameworkCoords, frameworkIsLocals, frameworkTags);
    }

    /***
    * MBX simulation box settings.
    ***/
    std::vector<double> box(9, 0.0);

    // Box x vector, used to be m11, m12, m13
    box[0] = simulationBox.cell.mm[0][0];
    box[1] = simulationBox.cell.mm[0][1];
    box[2] = simulationBox.cell.mm[0][2];

    // Box y vector, used to be m21, m22, m23
    box[3] = simulationBox.cell.mm[1][0];
    box[4] = simulationBox.cell.mm[1][1];
    box[5] = simulationBox.cell.mm[1][2];

    // Box z vector, used to be m31, m32, m33
    box[6] = simulationBox.cell.mm[2][0];
    box[7] = simulationBox.cell.mm[2][1];
    box[8] = simulationBox.cell.mm[2][2];
    
    for (int i = 0; i < 9; i++) {
        if (std::abs(box[i]) < 1e-10) {
            box[i] = 0.0; // Avoid numerical issues with very small values
        }
    }

    mbx->SetPBC(box);
    size_t numMon = mbx->GetNumMon();
    if (numMon <= 0)
    {
        RunningEnergy energySystem{};
        delete mbx;
        return energySystem;
    }
    /***
     * Calculate MBX energy.
    ***/

    bool do_grads = false;

    double e1b = mbx->OneBodyEnergy(do_grads);
    double e2b = mbx->TwoBodyEnergy(do_grads);
    double e3b = mbx->ThreeBodyEnergy(do_grads);
    double e4b = mbx->FourBodyEnergy(do_grads);
    double edisp = mbx->Dispersion(do_grads);
    double ebuck = mbx->Buckingham(do_grads);
    double elj = mbx->LennardJones(do_grads);
    mbx->Electrostatics(do_grads);
    double eelec_perm = mbx->GetPermanentElectrostaticEnergy();
    double eelec_ind = mbx->GetInducedElectrostaticEnergy();

    
    // double energy = e1b + e2b + e3b + e4b + edisp + ebuck + elj + eelec_perm + eelec_ind;    // full version
    double energy = e2b + e3b + e4b + edisp + eelec_perm + eelec_ind - system.elecPermFrameworkMBX; 
    // Excluding 1-body energy and the intra-molecular electrostatic energy of the framework

    // std::cerr << "e1b (omitted): " << e1b << std::endl;
    // std::cerr << "e2b: " << e2b << std::endl;
    // std::cerr << "e3b: " << e3b << std::endl;
    // std::cerr << "e4b: " << e4b << std::endl;
    // std::cerr << "edisp: " << edisp << std::endl;
    // // std::cerr << "ebuck: " << ebuck << std::endl;
    // // std::cerr << "elj: " << elj << std::endl;
    // std::cerr << "eelec_perm: " << eelec_perm << std::endl;
    // std::cerr << "eelec_ind: " << eelec_ind << std::endl;
    // std::cerr << "total_e: " << energy << std::endl;

    /***
     * Destroy MBX System object.
    ***/

    delete mbx;

    // Conversion kcal/mol -> RASPA internal units for energy, using Units module
    energy /= Units::EnergyToKCalPerMol;

    RunningEnergy energySystem{};

    // Assigning it to mbxEnergy
    energySystem.mbxEnergy = energy; 
    return energySystem;
}


// Use to calculate energies for the given fixed positions of framework, adsorbate molecules, and selected trial molecule for the MC move.
// The energy of the system can be calculated by either including or excluding the trial molecule.
RunningEnergy Interactions::computeMBXEnergy(
        const System &system,
        const std::vector<Component> &components,
        const SimulationBox &simulationBox,
        const std::optional<Framework> &framework,
        size_t selectedComponent,
        std::span<const Atom> frameworkAtoms,
        std::span<const Atom> moleculeAtoms,
        std::span<const Atom> selectedMoleculeAtoms,
        bool includeSelectedMoleculeAtoms,
        std::vector<double>* mbxEnergyLog
) noexcept
{

    /***
    * Setup MBX System object.
    ***/
    bblock::System* mbx = new bblock::System();
    int cumulativeTagIndex = 1;
    
    /***
    * Adding molecule monomers.
    ***/
    // Adding monomers molecules for all the components except the selectedComponents's selectedMolecule
    int numAtoms{0};
    for (size_t i = 0; i < components.size(); ++i)
    {
        numAtoms = static_cast<int>(components[i].atoms.size());           // Not sure whether MBX can use size_t
        std::vector<double> atomCoordinates(3 * numAtoms, 0.0);
        std::vector<std::string> atomNames(numAtoms, "none");
        std::string compName = components[i].name;

        // Index of first molecule of the given component.
        size_t idx_start = system.indexOfFirstMolecule(i);
        idx_start -= frameworkAtoms.size();                                 // Remove number of framework atoms
        bool keepMolecule{false};
        
        // Looping over all the molecules of given components
        for (size_t j = 0; j < system.numberOfMoleculesPerComponent[i]; ++j)
        {   
            for (const Atom &atom : selectedMoleculeAtoms)
            {
                if ( (j == static_cast<size_t>(atom.moleculeId)) && (selectedComponent == static_cast<size_t>(atom.componentId)) ) 
                {
                    keepMolecule = false;
                    break;
                }
                else { keepMolecule = true;}
            }
            if ( keepMolecule )
            {
                for(size_t k = 0; k < numAtoms; ++k) 
                {
                    atomCoordinates[k * 3] = moleculeAtoms[idx_start + k + j*numAtoms].position.x;
                    atomCoordinates[k * 3 + 1] = moleculeAtoms[idx_start + k + j*numAtoms].position.y;
                    atomCoordinates[k * 3 + 2] = moleculeAtoms[idx_start + k + j*numAtoms].position.z;
                    atomNames[k] = std::to_string(moleculeAtoms[idx_start + k + j*numAtoms].type);
                }

                // This decides weather this guest molecule is a real molecule or a ghost molecule
                // (such as an image of a molecule from a neighboring subdomain)
                // Since we use only one subdomain, we set this to true
                size_t islocal = 1;

                size_t tagIndex = cumulativeTagIndex;
                
                mbx->AddMonomer(atomCoordinates, atomNames, compName, islocal, tagIndex);

                cumulativeTagIndex += numAtoms;
            }
        }
    }
    if (includeSelectedMoleculeAtoms)
    {
        // Adding monomers for the selected component's selectedMolecule
        numAtoms = static_cast<int>(components[selectedComponent].atoms.size());    // Not sure whether MBX can use size_t
        std::vector<double> atomCoordinates(3 * numAtoms, 0.0);
        std::vector<std::string> atomNames(numAtoms, "none");
        std::string compName = components[selectedComponent].name;
        
        for(size_t k = 0; k < numAtoms; ++k) 
        {
            atomCoordinates[k * 3] = selectedMoleculeAtoms[k].position.x;
            atomCoordinates[k * 3 + 1] = selectedMoleculeAtoms[k].position.y;
            atomCoordinates[k * 3 + 2] = selectedMoleculeAtoms[k].position.z;
            atomNames[k] = std::to_string(selectedMoleculeAtoms[k].type);
        }
        size_t islocal = 1;
        size_t tagIndex = cumulativeTagIndex;
        mbx->AddMonomer(atomCoordinates, atomNames, compName, islocal, tagIndex);
        cumulativeTagIndex += numAtoms;
    }
    
    /***
    * MBX Miscellaneous settings.
    ***/
    // We now pass MBX file from RASPA.
    // char* json_path = "mbx.json";
    std::string json_path = system.mbxSettingsFilePath;
    std::ifstream t(json_path);
    t.seekg(0, std::ios::end);
    int size = t.tellg();
    std::string json_settings;
    json_settings.resize(size);
    t.seekg(0);
    t.read(&json_settings[0], size);
    mbx->SetUpFromJson(json_settings);

    
    /***
    * Adding framework atoms.
    ***/
    // Framework  --> In RASPA3 sometimes the framework is std::nullopt
    // Need to know MBX behaviour without the framework
    double eelec_perm_framework = 0.0;
    if (framework.has_value())
    {
        // eelec_perm_framework = Interactions::computeFrameworkElecPermMBXEnergy(system, simulationBox, framework, frameworkAtoms);

        std::vector<double> frameworkCoords(frameworkAtoms.size() * 3);
        std::vector<double> frameworkCharges(frameworkAtoms.size());
        std::vector<size_t> frameworkIsLocals(frameworkAtoms.size(), 1);
        std::vector<int> frameworkTags(frameworkAtoms.size());

        for (int i = 0; i < frameworkAtoms.size(); i++) {
            const Atom &atom = frameworkAtoms[i];

            frameworkCoords[i * 3] = atom.position.x;
            frameworkCoords[i * 3 + 1] = atom.position.y;
            frameworkCoords[i * 3 + 2] = atom.position.z;

            frameworkCharges[i] = atom.charge;              // TODO: MBX expects charges in units of e, but RASPA may give other units.
            frameworkTags[i] = cumulativeTagIndex;

            cumulativeTagIndex += 1;
        }
    
        mbx->SetExternalChargesAndPositions(frameworkCharges, frameworkCoords, frameworkIsLocals, frameworkTags);

    }

    /***
    * MBX simulation box settings.
    ***/
    std::vector<double> box(9, 0.0);

    // Box x vector, used to be m11, m12, m13
    box[0] = simulationBox.cell.mm[0][0];
    box[1] = simulationBox.cell.mm[0][1];
    box[2] = simulationBox.cell.mm[0][2];

    // Box y vector, used to be m21, m22, m23
    box[3] = simulationBox.cell.mm[1][0];
    box[4] = simulationBox.cell.mm[1][1];
    box[5] = simulationBox.cell.mm[1][2];

    // Box z vector, used to be m31, m32, m33
    box[6] = simulationBox.cell.mm[2][0];
    box[7] = simulationBox.cell.mm[2][1];
    box[8] = simulationBox.cell.mm[2][2];
    
    for (int i = 0; i < 9; i++) {
        if (std::abs(box[i]) < 1e-10) {
            box[i] = 0.0; // Avoid numerical issues with very small values
        }
    }

    mbx->SetPBC(box);
    size_t numMon = mbx->GetNumMon();
    if (numMon <= 0)
    {
        RunningEnergy energySystem{};
        delete mbx;
        return energySystem;
    }
    /***
     * Calculate MBX energy.
    ***/
    bool do_grads = false;

    double e1b = mbx->OneBodyEnergy(do_grads);
    double e2b = mbx->TwoBodyEnergy(do_grads);
    double e3b = mbx->ThreeBodyEnergy(do_grads);
    double e4b = mbx->FourBodyEnergy(do_grads);
    double edisp = mbx->Dispersion(do_grads);
    double ebuck = mbx->Buckingham(do_grads);
    double elj = mbx->LennardJones(do_grads);
    mbx->Electrostatics(do_grads);
    double eelec_perm = mbx->GetPermanentElectrostaticEnergy();
    double eelec_ind = mbx->GetInducedElectrostaticEnergy();

    // double energy = e1b + e2b + e3b + e4b + edisp + ebuck + elj + eelec_perm + eelec_ind;    // full version
    double energy = e2b + e3b + e4b + edisp + eelec_perm + eelec_ind - system.elecPermFrameworkMBX; 
    // Excluding 1-body energy and the intra-molecular electrostatic energy of the framework

    if (mbxEnergyLog)
    {
        if (mbxEnergyLog->size() != 7)
        {
            std::cerr << "mbxEnergyLog size mismatch" << "\n";
        }
        else{
            (*mbxEnergyLog)[0] = e1b;
            (*mbxEnergyLog)[1] = e2b;
            (*mbxEnergyLog)[2] = e3b;
            (*mbxEnergyLog)[3] = e4b;
            (*mbxEnergyLog)[4] = edisp;
            (*mbxEnergyLog)[5] = eelec_perm;
            (*mbxEnergyLog)[6] = eelec_ind;
        }
    }
    // std::cerr << "e1b (omitted): " << e1b << std::endl;
    // std::cerr << "e2b: " << e2b << std::endl;
    // std::cerr << "e3b: " << e3b << std::endl;
    // std::cerr << "e4b: " << e4b << std::endl;
    // std::cerr << "edisp: " << edisp << std::endl;
    // // std::cerr << "ebuck: " << ebuck << std::endl;
    // // std::cerr << "elj: " << elj << std::endl;
    // std::cerr << "eelec_perm: " << eelec_perm << std::endl;
    // std::cerr << "eelec_ind: " << eelec_ind << std::endl;
    // std::cerr << "total_e: " << energy << std::endl;

    // double energy =  mbx->Energy(false); // false means don't calculate forces
    // std::cerr << "MBX energy: " << energy << "[kcal/mol]" << std::endl;
    /***
     * Destroy MBX System object.
    ***/

    delete mbx;

    // Conversion kcal/mol -> RASPA internal units for energy, using Units module
    energy /= Units::EnergyToKCalPerMol;

    RunningEnergy energySum{};

    // Assigning it to mbxEnergy
    energySum.mbxEnergy = energy; 
    return energySum;
}

// used in translation.cpp, rotation.cpp
[[nodiscard]] RunningEnergy Interactions::computeMBXEnergyDifference(
        const System &system,
        const std::vector<Component> &components,
        const SimulationBox &simulationBox,
        const std::optional<Framework> &framework,
        size_t selectedComponent,
        std::span<const Atom> frameworkAtoms,
        std::span<const Atom> moleculeAtoms,
        std::span<const Atom> newatoms,
        std::span<const Atom> oldatoms
) noexcept
{        
    // Energy of the configuration before the trial MC move
    RunningEnergy oldEnergy = Interactions::computeMBXEnergy(system, components, simulationBox, framework, selectedComponent, frameworkAtoms, moleculeAtoms, oldatoms, true);
    
    // Energy of the configuration after the trial MC move
    RunningEnergy newEnergy = Interactions::computeMBXEnergy(system, components, simulationBox, framework, selectedComponent, frameworkAtoms, moleculeAtoms, newatoms, true);
    
    return newEnergy - oldEnergy;
}
