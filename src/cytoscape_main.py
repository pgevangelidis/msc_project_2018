# This program, reads all stored BGC objects from the src folder and reads the gamma vectors
import glob
import os.path
import numpy as np
import copy

if __name__ == '__main__':
    # The path to all BGC text files:
    path = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src\bgc_objects'
    # And this path contains the saved files for cytoscape
    save_path = r'C:\Users\user\Documents\msc_thesis_2018\msc_project_2018\src\cytoscape_files'

    output_style = "styleList_lda.csv"
    s_path = os.path.join(save_path, output_style)
    style_handler = open(s_path, 'w')
    ######
    # There are 3 different files for the cytoscape.
    # Each has a different threshold
    # Threshld A.(=1) B.(1 - 0.000005) C.(1 - 0.00001)
    output_file_a = 'bgc_topic_lda_1.csv'
    a_path = os.path.join(save_path, output_file_a)
    file_handler_a = open(a_path, 'w')

    output_file_b = 'bgc_topic_lda_999995.csv'
    b_path = os.path.join(save_path, output_file_b)
    file_handler_b = open(b_path, 'w')

    output_file_c = 'bgc_topic_lda_99999.csv'
    c_path = os.path.join(save_path, output_file_c)
    file_handler_c = open(c_path, 'w')

    file_handler_a.write('BGC,topic,threshold\n')
    file_handler_b.write('BGC,topic,threshold\n')
    file_handler_c.write('BGC,topic,threshold\n')
    style_handler.write('name,type\n')

    # loop through text files:
    numFiles = 0
    bgcName = ""
    print('Processing...\nPlease wait...\n')
    for filename in glob.glob(os.path.join(path, '*.txt')):
        bgc = open(filename, 'r')
        for line in bgc:
            if 'name:' in line:
                bgcName = line.split()[1]
                style_handler.write('{},BGC\n'.format(bgcName))
            if 'gamma:' in line:
                t = next(bgc).split()
                t = list(map(float, t))
                gamma = np.asarray(t)
                gamma = gamma/np.sum(gamma)
                for i in range(gamma.shape[0]):
                    if (gamma[i]>0.99999):
                        file_handler_c.write('{},topic{},{}\n'.format(bgcName,(i+1),gamma[i]))
                        if(gamma[i]>0.999995):
                            file_handler_b.write('{},topic{},{}\n'.format(bgcName,(i+1),gamma[i]))
                            if(gamma[i]==1):
                                file_handler_a.write('{},topic{},{}\n'.format(bgcName,(i+1),gamma[i]))
                break
        numFiles+=1
    for i in range(50):
        style_handler.write('topic{},topic\n'.format((i+1)))

    style_handler.close()
    file_handler_a.close()
    file_handler_b.close()
    file_handler_c.close()
    print('process done!')
