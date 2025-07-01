import os, shutil
import numpy as np
import argparse

def cleanup_directory(dir):
    for files in os.listdir(dir):
        path = os.path.join(dir, files)
        try:
            shutil.rmtree(path)
        except OSError:
            os.remove(path)


def modify_imaging_script_list(imaging_script_list, imagedir, sequence, temperature, frame_ind, movie=True):
    if movie:
        modified_imaging_script_list = [s.replace('GREEN_FILENAME', get_image_filename_movie(imagedir, sequence, temperature, 'Green', frame_ind)).replace('RED_FILENAME', get_image_filename_movie(imagedir, sequence, temperature, 'Red', frame_ind)) \
            for s in imaging_script_list]
    else:
        modified_imaging_script_list = [s.replace('GREEN_FILENAME', get_image_filename(imagedir, sequence, temperature, 'Green')).replace('RED_FILENAME', get_image_filename(imagedir, sequence, temperature, 'Red')) \
            for s in imaging_script_list]

    return modified_imaging_script_list

# ====== functions to get image filename ======
def get_image_filename_movie(imagedir, sequence, temperature, channel, frame_ind):
    """
    channel - 'Green' or 'Red'
    """
    filename = '\"%s\\%s%02d-%d_%d\\NNNlib2b_DNA_tileVAR_tile_%s_600ms_VAR_TIMESTAMP.tif\"' % (imagedir, channel, sequence, frame_ind, temperature, channel.lower())
    return filename

def get_image_filename(imagedir, sequence, temperature, channel):
    """
    channel - 'Green' or 'Red'
    """
    filename = '\"%s\\%s%02d_%.1f\\NNNlib2b_DNA_tileVAR_tile_%s_600ms_VAR_TIMESTAMP.tif\"' % (imagedir, channel, sequence, temperature, channel.lower())
    return filename

# ====== functions to make folders =====

def make_folders_movie(folderdir, seq, frame_ind, temperature):
    os.makedirs(os.path.join(folderdir, '%02d-%d_Green%d' % (seq, frame_ind, temperature)))
    os.makedirs(os.path.join(folderdir, '%02d-%d_Red%d' % (seq, frame_ind, temperature)))

def make_folders(folderdir, seq, temperature):
    os.makedirs(os.path.join(folderdir, 'Green%02d_%s' % (seq, '{:g}'.format(temperature))))
    os.makedirs(os.path.join(folderdir, 'Red%02d_%s' % (seq, '{:g}'.format(temperature))))

# ====== make scripts and folders, calls functions above ======
def make_scripts_folders(first_melt_sequence, temperatures, n_frame, array_image_dir, template_file, script_outdir, folder_outdir):
    
    n_temperature = len(temperatures)
    sequences = np.arange(first_melt_sequence, first_melt_sequence + n_temperature)
    experiment_dict = dict(zip(sequences, temperatures))

    cleanup_directory(folder_outdir)
    cleanup_directory(script_outdir)

    with open(template_file, 'r') as fh:
        template_content = fh.readlines()

    start_loop_ind = template_content.index('<sub name="imagetiles">\n')
    end_loop_ind = template_content.index('</sub>\n')
    imaging_script_list = template_content[start_loop_ind:end_loop_ind+1]

    for seq, t in experiment_dict.items():
        
        if n_frame > 1:
            movie_imaging_script_list = []
            for frame_ind in range(1, n_frame + 1):
                movie_imaging_script_list = movie_imaging_script_list + modify_imaging_script_list(imaging_script_list, array_image_dir, sequence=seq, temperature=t, frame_ind=frame_ind) 
                make_folders_movie(folder_outdir, seq, frame_ind, t)
        elif n_frame == 1:
            movie_imaging_script_list = modify_imaging_script_list(imaging_script_list, array_image_dir, sequence=seq, temperature=t, frame_ind=1, movie=False)
            make_folders(folder_outdir, seq, t)

        outfile = os.path.join(script_outdir, '%02d_ImageBothTiles_%.1f.xml' % (seq, t))
        with open(outfile, 'w+') as fh:
            fh.writelines(template_content[:start_loop_ind])
            fh.writelines(movie_imaging_script_list)
            fh.writelines([s.replace('TEMPERATURE', '%.1f'%t) for s in template_content[end_loop_ind + 1:]])

if __name__ == "__main__":
    # I'm too lazy to write boiler plate argparse
    # parser = argparse.ArgumentParser()
    # parser.add_argument('-fs', '--first_melt_sequence', type=int, help='first step number when the melt curve starts', required=True)
    # parser.add_argument('-ft')

    # args = parser.parse_args()


    # ====== Settings =======
    # First step number when the melt curve starts
    first_melt_sequence = 1
    # Define the timerature points
    temperatures = np.arange(20, 62.5, 2.5)
    temperatures = np.append(temperatures, 20)
    # Number of frames to take, usually 1 
    n_frame = 1
    # The data folder on the array computer, used for changing filenames in the xml script
    array_image_dir = 'D:\\images\\NNN-melt\\NNN_20211222'

    # The template script where variables are changed to GREEN_FILENAME, TEMPERATURE, etc.
    template_file = r'../data/imaging_scripts/ImageBothTiles_TEMPLATE.xml'
    # Where to put the output imaging scripts
    script_outdir = r'../data/NNN_20211222/imaging_scripts/'
    # Path to put the empty folders
    folder_outdir = r'../data/NNN_20211222/imaging_folders/'
    # =======================
    
    make_scripts_folders(first_melt_sequence, temperatures, n_frame, array_image_dir, template_file, script_outdir, folder_outdir)