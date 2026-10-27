"""`Rscript` on the agent PATH: run the host R with the environment it needs.

The acumen sandbox blanks every variable outside its allowlist, so the shim sets
R_HOME, the R shared-library path and the user library itself. The defaults are
this host's R 4.5.2 module; override each with LUCID_R_HOME / LUCID_R_LIBS_USER.
"""

import os
import sys

# ponytail: host paths baked in; use env_passthrough once a second host runs it.
R_HOME = os.environ.get(
    "LUCID_R_HOME", "/share/software/statistical/R/R-4.5.2/gcc-8.5.0/R-4.5.2/lib64/R"
)
R_LIBS_USER = os.environ.get(
    "LUCID_R_LIBS_USER", "/exports/para-lipg-hpc/mdmanurung/R/4.5"
)
# Newer libstdc++ than the system's: user-library farver needs GLIBCXX_3.4.29.
GCC_RUNTIME = os.environ.get(
    "LUCID_GCC_RUNTIME", "/home/mdmanurung/.local/lib/conda-gcc-runtime"
)


def main():
    env = dict(os.environ, R_HOME=R_HOME, R_LIBS_USER=R_LIBS_USER)
    env["LD_LIBRARY_PATH"] = os.pathsep.join(
        filter(None, [f"{R_HOME}/lib", GCC_RUNTIME, env.get("LD_LIBRARY_PATH")])
    )
    rscript = os.path.join(R_HOME, "bin", "Rscript")
    os.execve(rscript, [rscript, *sys.argv[1:]], env)
