# Reproducing the bounded checks

The local candidate source is `development/RASPA3-MBX-r153-tests`, branch `review/r153-general-tests`. The base is published df7a9065. `repairs.patch` contains the exact six-file change and is suitable for review or applying to that base in another isolated checkout. Do not apply it to an active calculation/deployment directory.

`configure_final_command.json` and `configure_off_final_command.json` contain the successful CMake configuration arguments for this machine; use the recorded Clang20.1.6/conda/MBX dependency paths. They include CONDA=TRUE and LLD, omitted in the preserved initial configuration attempt. To run one configuration without shell interpolation:

```python
import json, os, subprocess
from pathlib import Path
report = Path('/home/yuc/calculation/RASPA-MBX/code_review/20261023_r153')
tc = '/opt/ohpc/pub/apps/miniconda3/envs/raspa3-mbx'
env = dict(os.environ, PATH=tc+'/bin:'+os.environ['PATH'], LD_LIBRARY_PATH=tc+'/lib')
cmd = json.loads((report/'configure_final_command.json').read_text())
subprocess.run(cmd, env=env, check=True)
subprocess.run([tc+'/bin/cmake','--build',cmd[cmd.index('-B')+1],
                '--target','unit_tests_raspakit','unit_tests_mbx','-j4'], env=env, check=True)
```

The MBX-OFF target list omits unit_tests_mbx. The complete build logs, caches and runtime identities are retained. From the project root, native checks use the following commands; the final label creates a separate output directory and prevents overwriting previous receipts:

```sh
python3 code_review/20261023_r153/run_validation.py on full repeat1
python3 code_review/20261023_r153/run_validation.py on mbx repeat1
python3 code_review/20261023_r153/run_validation.py off affected repeat1
```

Each command runs in its own NFS output directory with one OpenMP/OpenBLAS thread, writes log/XML/JSON and exits nonzero on failure. `baseline six` runs the unchanged original published executable. Do not treat CTest source-scanned registrations as actual executed cases. Standalone diagnostic sources and commands live in hessian/, hybrid_io/, vdw/ and sphere/. They do not write original simulation data or submit jobs.
