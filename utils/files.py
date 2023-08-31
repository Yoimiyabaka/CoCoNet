import os
import shutil


def move_file_to_dir(file_path, target_dir, copy=True):
    # Determine the file name
    file_name = os.path.basename(file_path)

    # Determine the destination path
    destination_path = os.path.join(target_dir, file_name)
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # Move the file to the target directory
    if copy:
        shutil.copy2(file_path, destination_path)
    else:
        shutil.move(file_path, destination_path)


def move_dir_to_dir(source_dir, target_dir, copy=True):
    '''
        move all files under one directory to another
    '''
    # Create the destination folder if it doesn't exist
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # Get a list of all files in the source folder
    files = os.listdir(source_dir)

    # Move each file to the destination folder
    for file in files:
        source_path = os.path.join(source_dir, file)
        target_path = os.path.join(target_dir, file)
        if copy:
            shutil.copy2(source_path, target_path)
        else:
            shutil.move(source_path, target_path)


def validate_input(path):
    if os.path.isdir(path):
        if is_folder_with_valid_fasta_files(path):
            return "Folder"
    elif os.path.isfile(path) and path.lower().endswith('.fasta'):
        return "FASTA file"
    else:
        raise ValueError("Invalid input: " + path)


def is_folder_with_valid_fasta_files(folder_path):

    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        if os.path.isfile(file_path):
            if file_name.lower().endswith('.fasta'):
                # Validate .fasta file
                continue
            elif file_name.lower().endswith('.fasta.gz'):
                # Validate .fasta.gz file
                continue
            else:
                # Invalid file extension
                raise ValueError("Invalid FASTA file: " + file_name)
    return True

