import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cv2
import os

def Visualization(Tinit, TField, info):
    
    folder = 'figures_'+Tinit+'/'
    i = 0
    print('Creating figures for each timestep '+Tinit+' start (.png) ...')
    for j in np.arange(0, len(Tfield)-info.ny[0], info.ny[0]):
        plt.figure(figsize=((32, 8)))
        plt.contourf(Tfield[j:j+info.ny[0], :], np.linspace(0, 1, 11))
        plt.title('Temperature field - '+Tinit+' start - Pr =%f' %info['Pr'])    
        plt.gca().set_aspect("equal")
        plt.colorbar(shrink=0.2)
        if i < 10 :
            plt.savefig(folder + 'TField000'+str(i)+'.png',bbox_inches='tight')
        elif i < 100 :
            plt.savefig(folder + 'TField00'+str(i)+'.png',bbox_inches='tight')
        elif i < 1000 :
            plt.savefig(folder + 'TField0'+str(i)+'.png',bbox_inches='tight')
        else :
            plt.savefig(folder + 'TField'+str(i)+'.png',bbox_inches='tight')
        plt.clf()
        plt.close()
        i += 1

    print('Figures successfully created !!')
    print('Generating a movie for the '+Tinit+' start ...')
    image_folder = 'figures_'+Tinit
    video_name = 'convection_'+Tinit+'_Pr'+str(round(float(info.Pr), 2))+'.avi'
    images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, layers = frame.shape

    video = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'DIVX'), 8, (width,height))

    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))

    video.release()
    print('Movie successfully created !!')

# Define input and call function
Tinit = 'random'
Tfield = np.loadtxt('results/Tfield_'+Tinit+'.txt')
info = pd.read_csv('results/info_'+Tinit+'.csv')
Visualization(Tinit, Tfield, info)