file(MAKE_DIRECTORY "${WORK}")
file(COPY "${FIXTURE}/co2.json" "${FIXTURE}/force_field.json" "${FIXTURE}/mbx.json" DESTINATION "${WORK}")
set(input [=[{
  "SimulationType": "MonteCarlo", "NumberOfProductionCycles": 1, "NumberOfInitializationCycles": 0,
  "NumberOfEquilibrationCycles": 0, "WriteBinaryRestartEvery": 1, "RandomSeed": 42, "ForceField": ".",
  "Systems": [{"Type": "Box", "BoxLengths": [28,28,28], "ExternalTemperature": 298.0,
               "UseMBX": true, "MBXSettingsFile": "mbx.json"}],
  "Components": [{"Name": "co2", "MoleculeDefinition": ".", "TranslationProbability": 1.0,
                  "CreateNumberOfMolecules": 1}]
}]=])
file(WRITE "${WORK}/simulation.json" "${input}")
execute_process(COMMAND "${RASPA}" WORKING_DIRECTORY "${WORK}" RESULT_VARIABLE status
                OUTPUT_QUIET ERROR_VARIABLE diagnostic)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Could not create the MBX checkpoint: ${diagnostic}")
endif()
string(JSON input SET "${input}" RestartFromBinaryFile true)
string(JSON input SET "${input}" Systems 0 UseMBX false)
file(WRITE "${WORK}/simulation.json" "${input}")
execute_process(COMMAND "${RASPA}" WORKING_DIRECTORY "${WORK}" RESULT_VARIABLE status
                OUTPUT_QUIET ERROR_VARIABLE diagnostic)
if(status EQUAL 0 OR NOT diagnostic MATCHES "does not match the energy model")
  message(FATAL_ERROR "Expected a model mismatch when reading the MBX checkpoint: ${diagnostic}")
endif()
