# Computer setup

[Back to the three-step guide](README.md)

Use a 64-bit Linux computer with an Intel or AMD processor. The commands below
are for **Ubuntu 22.04 or newer**. If you use a shared university server, ask
its administrator whether Python, Git, and Conda are already available.

## Install Python and Git

Open Terminal and run:

```bash
sudo apt update
sudo apt install -y python3 git curl
```

Your computer may ask for your login password. Check the installation:

```bash
python3 --version
git --version
```

Python must be version 3.10 or newer.

## Install Conda

If you already have Conda, Miniforge, Miniconda, or Anaconda installed, skip this
step. Otherwise, download and install Miniforge:

```bash
curl -fL -o Miniforge3.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
```

The reproduction program can find this installation automatically. You do not
need to activate an environment or change your terminal settings.

These commands use the project's [official Linux installation instructions](https://github.com/conda-forge/miniforge#unix-like-platforms-macos-linux--wsl).

## Choose a place for the data

Use a Linux drive with plenty of free space. Plan for 700 GB in the project
folder's drive and 100 GB in `/tmp`, or 800 GB when both share one drive. Keep
at least 32 GiB RAM available; a machine with 64 GB total RAM gives more room.
These are planning allowances, not measured minimum requirements.

The program checks the actual locations and tells you if something is missing.
Return to [step 2 of the main guide](README.md#2-download-this-repository).
