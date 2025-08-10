import os
def setup_directories(path):
    os.makedirs(f"{path}", exist_ok=True)
    os.makedirs(f"{path}/img", exist_ok=True)
    os.makedirs(f"{path}/data", exist_ok=True)