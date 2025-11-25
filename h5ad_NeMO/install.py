import subprocess
import sys
import importlib

REQUIRED_PACKAGES = [
    'scipy',
    'tqdm',
    'openpyxl',
    'h5py',
    'sklearn',
    'scanpy',
    'pandas',
]

def check_package(package_name):
    """Check if package is already installed"""
    try:
        importlib.import_module(package_name)
        return True
    except ImportError:
        return False

def install_package(package_name):
    """Install single package"""
    print(f"Installing {package_name}...")
    try:
        subprocess.check_call([
            sys.executable, "-m", "pip", "install", "--user", package_name
        ], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        print(f"{package_name} installed successfully")
        return True
    except Exception as e:
        print(f"\033[91m{package_name}\033[0m installation failed: {e}")
        return False

def main():
    print("=" * 50)
    print("Checking and installing program dependencies")
    print("=" * 50)
    to_install = []
    for package in REQUIRED_PACKAGES:
        if check_package(package):
            print(f"{package} already installed")
        else:
            installed = install_package(package)
            if not installed:
                to_install.append(package)
    if len(to_install) > 0:
        print("The following dependencies failed to install, please install manually")
        print(f"\033[91m{' '.join(to_install)}\033[0m")
    else:
        print("All dependencies installed successfully!")

if __name__ == "__main__":
    main()
