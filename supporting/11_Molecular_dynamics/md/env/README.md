# Software environment for the molecular-dynamics set-up

Prefix on the Pi: `/root/cqc-sfa-tools/md` (aarch64, conda-forge), made on 5 Oct 2026 with

    micromamba create -y -p /root/cqc-sfa-tools/md -c conda-forge python=3.12 ambertools=26 openmm mdanalysis parmed numpy scipy matplotlib

`md_env.yml` lists the packages with versions, `md_env_explicit.txt` the exact files (rebuild with
`micromamba create -p <prefix> --file md_env_explicit.txt`). The quantum-chemical steps use the PySCF environment of
the DFT part (`/root/cqc-sfa-tools/env`, PySCF 2.14.0), the pose analysis the docking environment
(`/root/cqc-sfa-tools/dock`). OpenMM offers only the Reference and CPU platforms on this machine.
