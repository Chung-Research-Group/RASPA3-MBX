module;

module system;

import std;
import atom;

bool System::mbxAvailable() { return false; }

void System::configureMBX(const std::string&)
{
  throw std::runtime_error("UseMBX requires RASPA built with -DBUILD_MBX=ON");
}

double System::mbxEnergy(std::size_t, std::optional<std::size_t>, std::span<const Atom>, std::array<double, 7>*) const
{
  throw std::runtime_error("This binary was built without MBX support");
}
