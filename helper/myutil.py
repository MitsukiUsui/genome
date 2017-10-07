import subprocess


def myrun(cmd):
    """
    run command using subprocess module
    """
    try:
        subprocess.run(cmd.split(), check=True)
    except subprocess.CalledProcessError as e:
        print(e)
        return False
    return True
