from openmm import Platform

# Try to get the fastest available platform
def get_best_platform():
    try:
        platform = Platform.getPlatformByName("CUDA")
        print("Using platform: CUDA (NVIDIA GPU)")
    except Exception:
        try:
            platform = Platform.getPlatformByName("OpenCL")
            print("Using platform: OpenCL")
        except Exception:
            platform = Platform.getPlatformByName("CPU")
            print("Using platform: CPU (no GPU found)")
    return platform