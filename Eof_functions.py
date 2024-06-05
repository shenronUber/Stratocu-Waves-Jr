from openpiv import tools, pyprocess, validation, filters, scaling 

import numpy as np
import cv2
import os
import glob

import matplotlib.pyplot as plt
from matplotlib import colormaps
from scipy.linalg import eigh
from eofs.standard import Eof
from scipy.fft import fft2, ifft2, fftn, ifftn, fftshift, ifftshift
from scipy.optimize import curve_fit

from PIL import Image
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter as smoo

from skimage import measure
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression

from mpl_toolkits.mplot3d import Axes3D

def Calculate_Advection(frame_a,frame_b):
    # High pass the images to make sure the tracking is not of wave velocity itself
    frame_a_hp = frame_a.copy() - smoo(frame_a, 10)
    frame_b_hp = frame_b.copy() - smoo(frame_b, 10)

    winsize = 36 # pixels, interrogation window size in frame A
    searchsize = 36  # pixels, search in image B big enough to contain credible velocity 
    overlap = 18 # pixels, 50% overlap if half of winsize
    dt = 1 # time interval between images, converts pixel displacement to velocity

    # Coordinates of velocity positions in image array
    x, y = pyprocess.get_coordinates(image_size=frame_a.shape, 
                                    search_area_size=searchsize, 
                                    overlap=overlap )

    # Use high-pass images as the PIV source frame_a_hp, frame_b_hp

    u, v, sig2noise = pyprocess.extended_search_area_piv(  frame_a_hp.astype(np.int32), 
                                                        frame_b_hp.astype(np.int32), 
                                                        window_size=winsize, 
                                                        overlap=overlap, 
                                                        dt=dt, 
                                                        search_area_size=searchsize, 
                                                        sig2noise_method='peak2mean')
    # return u/dt, v/dt, sig2noise

    mask = validation.global_std(u, v, std_threshold=3)
    #replace outliers with NaNs
    u[mask] = np.nan
    v[mask] = np.nan

    mask = validation.global_std(u, v, std_threshold=3)
    #replace outliers with NaNs
    u[mask] = np.nan
    v[mask] = np.nan

    u, v = filters.replace_outliers( u, v, mask,
                                    method='localmean', 
                                    max_iter=7, 
                                    kernel_size=7)
    
    x, y, u, v = tools.transform_coordinates(x, y, u, v)

    # sum/num will be the average, just remember to avoid zero division
    a2sum = frame_a.copy()*0.0
    a2num = frame_a.copy()*0

    # u and v offsets in nearest pixel, integer
    xi = x.round().astype(int).ravel()
    yi = y.round().astype(int).ravel()
    ui = u.round().astype(int).ravel()
    vi = v.round().astype(int).ravel()


    # Loop over all the little windows in a and add them into a2 
    for idx,vpoint in enumerate(xi.ravel()):
        windowa = (slice(yi[idx]-winsize//2, yi[idx]+winsize//2,1), 
                slice(xi[idx]-winsize//2, xi[idx]+winsize//2,1) ) # source image

        windowa2= (slice(yi[idx]-winsize//2 + -vi[idx], 
                        yi[idx]+winsize//2 + -vi[idx],1), 
                slice(xi[idx]-winsize//2 + ui[idx], 
                        xi[idx]+winsize//2 + ui[idx],1) ) # where to put it
        
        try: 
            a2sum[windowa2] += frame_a[windowa]
            a2num[windowa2] += 1
        except: 
            continue
            
    a2 = (a2sum/(a2num 
                +0.000001)) # avoid division by zero
    return a2

def Show_Advection(frame_a, frame_b):
    
    plt.figure(figsize=[10,5]) # type: ignore
    plt.subplot(1, 3, 1)
    plt.imshow(frame_a); plt.title('original frame'); plt.clim(0,200)
    plt.subplot(1, 3, 2)
    plt.imshow(Calculate_Advection(frame_a, frame_b)); plt.title('advected frame'); plt.clim(0,200)
    plt.subplot(1, 3, 3)
    plt.imshow(frame_b); plt.title('original next frame'); plt.clim(0,200)

def Estimate_Advection(frame_a,frame_b):

    plt.figure(figsize=[10,5]) # type: ignore
    plt.subplot(1, 2, 1)
    plt.imshow(Calculate_Advection(frame_a, frame_b)); plt.title('advected frame'); plt.clim(0,200); plt.colorbar(shrink=0.5)

    plt.subplot(1, 2, 2)
    plt.imshow(frame_b-Calculate_Advection(frame_a, frame_b),cmap='bwr'); plt.title('Next Frame - Advected Frame'); plt.clim(-50,50); plt.colorbar(shrink=0.5)

def distance_from(a, index):
    i,j = np.indices(a.shape, sparse=True)
    return np.sqrt((i-index[0])**2 + (j-index[1])**2)

def angle_from(a, index):
    i,j = np.indices(a.shape, sparse=True)
    return np.arctan2((j-index[1]),(i-index[0]))

def Show_Brightness(frame_a,frame_b):
    # Fourier transform and shift to center and total wavenumber array kl
    diff = frame_b-frame_a
    diffhat = np.fft.fftshift(np.fft.fft2(diff))

    # Power
    power = np.abs(diffhat)**2

    # Wavenumbers 
    kl = distance_from(diffhat, [diffhat.shape[0]/2, diffhat.shape[1]/2] )

    mask = (  (kl>0) & (kl<=10) ) 

    recon = np.fft.ifft2( np.fft.ifftshift( diffhat*mask ))

    plt.imshow(recon.real, cmap='RdBu_r')
    plt.clim(-20,20); plt.colorbar(shrink=0.7); plt.title('filtered d/dt(brightness)')

def Show_Convergence(frame_a,frame_b,SmoothBorders=False):
    # High pass the images to make sure the tracking is not of wave velocity itself
    frame_a_hp = frame_a.copy() - smoo(frame_a, 10)
    frame_b_hp = frame_b.copy() - smoo(frame_b, 10)

    winsize = 36 # pixels, interrogation window size in frame A
    searchsize = 36  # pixels, search in image B big enough to contain credible velocity 
    overlap = 18 # pixels, 50% overlap if half of winsize
    dt = 1 # time interval between images, converts pixel displacement to velocity

    # Coordinates of velocity positions in image array
    x, y = pyprocess.get_coordinates(image_size=frame_a.shape, 
                                    search_area_size=searchsize, 
                                    overlap=overlap )

    # Use high-pass images as the PIV source frame_a_hp, frame_b_hp

    u, v, sig2noise = pyprocess.extended_search_area_piv(  frame_a_hp.astype(np.int32), 
                                                        frame_b_hp.astype(np.int32), 
                                                        window_size=winsize, 
                                                        overlap=overlap, 
                                                        dt=dt, 
                                                        search_area_size=searchsize, 
                                                        sig2noise_method='peak2mean')

    mask = validation.global_std(u, v, std_threshold=3)
    #replace outliers with NaNs
    u[mask] = np.nan
    v[mask] = np.nan

    mask = validation.global_std(u, v, std_threshold=3)
    #replace outliers with NaNs
    u[mask] = np.nan
    v[mask] = np.nan

    u, v = filters.replace_outliers( u, v, mask,
                                    method='localmean', 
                                    max_iter=7, 
                                    kernel_size=7)
    
    x, y, u, v = tools.transform_coordinates(x, y, u, v)

    divx = np.gradient(u)[1]
    divy = -np.gradient(v)[0]
    div = divx+divy
    if SmoothBorders == True:
        window = np.hanning(8)  # 8 pixels wide hanning window (4*2)

        # Apply the window to top and bottom borders
        div[:4, :] = div[:4, :] * window[:4, None]
        div[-4:, :] = div[-4:, :] * window[:4, None]

        # Apply the window to left and right borders
        div[:, :4] = div[:, :4] * window[:4, None].T
        div[:, -4:] = div[:, -4:] * window[:4, None].T

        plt.imshow(smoo(-div, 2), cmap='RdBu_r', clim=[-0.2,0.2]); plt.colorbar(shrink=0.7)
        plt.title('convergence of PIV velocities')  
    else:
        plt.imshow(smoo(-div, 2), cmap='RdBu_r', clim=[-1,1]); plt.colorbar(shrink=0.7)
        plt.title('convergence of PIV velocities')    

def Show_Contour(frame_a,frame_b,SmoothBorders=False):
    # High pass the images to make sure the tracking is not of wave velocity itself
    frame_a_hp = frame_a.copy() - smoo(frame_a, 10)
    frame_b_hp = frame_b.copy() - smoo(frame_b, 10)

    winsize = 36 # pixels, interrogation window size in frame A
    searchsize = 36  # pixels, search in image B big enough to contain credible velocity 
    overlap = 18 # pixels, 50% overlap if half of winsize
    dt = 1 # time interval between images, converts pixel displacement to velocity

    # Use high-pass images as the PIV source frame_a_hp, frame_b_hp

    u, v, sig2noise = pyprocess.extended_search_area_piv(  frame_a_hp.astype(np.int32), 
                                                        frame_b_hp.astype(np.int32), 
                                                        window_size=winsize, 
                                                        overlap=overlap, 
                                                        dt=dt, 
                                                        search_area_size=searchsize, 
                                                        sig2noise_method='peak2mean')
    # Coordinates of velocity positions in image array
    x, y = pyprocess.get_coordinates(image_size=frame_a.shape, 
                                    search_area_size=searchsize, 
                                    overlap=overlap )

    mask = validation.global_std(u, v, std_threshold=3)
    #replace outliers with NaNs
    u[mask] = np.nan
    v[mask] = np.nan

    mask = validation.global_std(u, v, std_threshold=3)
    #replace outliers with NaNs
    u[mask] = np.nan
    v[mask] = np.nan

    u, v = filters.replace_outliers( u, v, mask,
                                    method='localmean', 
                                    max_iter=7, 
                                    kernel_size=7)
    
    x, y, u, v = tools.transform_coordinates(x, y, u, v)

    divx = np.gradient(u)[1]
    divy = -np.gradient(v)[0]
    div = divx+divy
    if SmoothBorders == True:
        # transform the first 4 rows and columuns into NaN
        window = np.hanning(8)  # 8 pixels wide hanning window (4*2)

        # Apply the window to top and bottom borders
        div[:4, :] = div[:4, :] * window[:4, None]
        div[-4:, :] = div[-4:, :] * window[:4, None]

        # Apply the window to left and right borders
        div[:, :4] = div[:, :4] * window[:4, None].T
        div[:, -4:] = div[:, -4:] * window[:4, None].T


    # Fourier transform and shift to center and total wavenumber array kl
    diff = frame_b-frame_a
    diffhat = np.fft.fftshift(np.fft.fft2(diff))

    # Wavenumbers 
    kl = distance_from(diffhat, [diffhat.shape[0]/2, diffhat.shape[1]/2] )

    mask = (  (kl>0) & (kl<=10) ) 

    recon = np.fft.ifft2( np.fft.ifftshift( diffhat*mask ))

    plt.pcolormesh(recon.real, cmap='RdBu_r')
    plt.clim(-20,20); plt.title('d/dt(brightness) and conv(PIV wind)')

    plt.contour(x,y, np.flipud(smoo(-div, 2)), levels=(-5+np.arange(11))/20., cmap='RdBu_r')
    if SmoothBorders == False:
        plt.ylim(950,50); plt.xlim(500,1400)  
    
def Get_Matrices(frame_a,frame_b):
    # High pass the images to make sure the tracking is not of wave velocity itself
    frame_a_hp = frame_a.copy() - smoo(frame_a, 10)
    frame_b_hp = frame_b.copy() - smoo(frame_b, 10)

    winsize = 36 # pixels, interrogation window size in frame A
    searchsize = 36  # pixels, search in image B big enough to contain credible velocity 
    overlap = 18 # pixels, 50% overlap if half of winsize
    dt = 1 # time interval between images, converts pixel displacement to velocity

    # Use high-pass images as the PIV source frame_a_hp, frame_b_hp

    u, v, sig2noise = pyprocess.extended_search_area_piv(  frame_a_hp.astype(np.int32), 
                                                        frame_b_hp.astype(np.int32), 
                                                        window_size=winsize, 
                                                        overlap=overlap, 
                                                        dt=dt, 
                                                        search_area_size=searchsize, 
                                                        sig2noise_method='peak2mean')
    # Coordinates of velocity positions in image array
    x, y = pyprocess.get_coordinates(image_size=frame_a.shape, 
                                    search_area_size=searchsize, 
                                    overlap=overlap )

    mask = validation.global_std(u, v, std_threshold=3)
    #replace outliers with NaNs
    u[mask] = np.nan
    v[mask] = np.nan

    mask = validation.global_std(u, v, std_threshold=3)
    #replace outliers with NaNs
    u[mask] = np.nan
    v[mask] = np.nan

    u, v = filters.replace_outliers( u, v, mask,
                                    method='localmean', 
                                    max_iter=7, 
                                    kernel_size=7)
    
    x, y, u, v = tools.transform_coordinates(x, y, u, v)

    divx = np.gradient(u)[1]
    divy = -np.gradient(v)[0]
    div = divx+divy

    # Fourier transform and shift to center and total wavenumber array kl
    diff = frame_b-frame_a
    diffhat = np.fft.fftshift(np.fft.fft2(diff))

    # Wavenumbers 
    kl = distance_from(diffhat, [diffhat.shape[0]/2, diffhat.shape[1]/2] )

    mask = (  (kl>0) & (kl<=10) ) 

    recon = np.fft.ifft2( np.fft.ifftshift( diffhat*mask ))

    return smoo(-div, 2), recon.real

def resize_matrix(original, target_dimensions):

    # Determine the current dimensions of the original matrix
    original_dimensions = original.shape

    # Step 1: Define the original grid
    x = np.linspace(0, 1, original_dimensions[1])  
    y = np.linspace(0, 1, original_dimensions[0])  

    # Step 2: Define the function on this grid
    X, Y = np.meshgrid(x, y, indexing='ij')
    Z = original.T

    # Step 3: Create a finer grid for interpolation
    x_fine = np.linspace(0, 1, target_dimensions[1])  
    y_fine = np.linspace(0, 1, target_dimensions[0])  

    # Step 4: Use RegularGridInterpolator
    interpolator = RegularGridInterpolator((x, y), Z)
    X_fine, Y_fine = np.meshgrid(x_fine, y_fine, indexing='ij')
    points_fine = np.array([X_fine.ravel(), Y_fine.ravel()]).T
    Z_fine = interpolator(points_fine).reshape(target_dimensions[1], target_dimensions[0])

    return Z_fine.T

def Show_Candidates(frame_a,frame_b,plot_option='no'):
    con, bri = Get_Matrices(frame_a,frame_b)
    con=resize_matrix(con, bri.shape)

    # create a tuple array, that takes for each cell the con and bri values
    candidates = np.zeros((con.shape[0], con.shape[1]), dtype=[('con', float), ('bri', float)])
    candidates['con'] = con
    candidates['bri'] = bri

    # Define a threshold for closeness to the 1:1 line
    threshold = 0.5  

    # Get the points close to the 1:1 line
    close_to_line = np.abs(candidates['bri'] - candidates['con']) < threshold

    if plot_option == 'yes':
        # Scatter plot of con vs bri with the required axis settings and aspect ratio
        plt.figure(figsize=(10, 5))  

        #plot the close to line in red and the rest in blue
        plt.scatter(candidates['bri'][close_to_line], candidates['con'][close_to_line], s=1, color='red')
        plt.scatter(candidates['bri'][~close_to_line], candidates['con'][~close_to_line], s=1, color='blue')
        plt.xlabel('Brightness')
        plt.ylabel('Convergence')

        # Adjust axes to be centered on (0, 0)
        plt.axhline(0, color='black', linewidth=0.8)
        plt.axvline(0, color='black', linewidth=0.8)

        # Set the axis limits 
        plt.xlim(-10, 10)
        plt.ylim(-5, 5)

        # Plot a line with a slope of 1:1 for convergence:brightness
        plt.plot([bri.min(), bri.max()], [bri.min(), bri.max()], 'k--', lw=2)  # Plot a 1:1 line for reference

        plt.title('Scatter Plot of Convergence vs. Brightness')
        plt.show()

    # Retrieve the coordinates of the points close to the 1:1 line
    close_points_indices = np.nonzero(close_to_line)
    close_points_coordinates = list(zip(close_points_indices[0], close_points_indices[1]))

    # plot frame_b with the close points in red
    plt.figure(figsize=(10, 5))
    plt.imshow(frame_b, cmap='gray')
    plt.scatter([p[1] for p in close_points_coordinates], [p[0] for p in close_points_coordinates], s=1, color='red')
    plt.title('Potential candidates for GW')
    plt.show()

def save_frames_from_video(video_path, output_folder):
    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    cap = cv2.VideoCapture(video_path)
    if not cap.isOpened():
        print("Error: Could not open video.")
        return

    count = 0
    while cap.isOpened():
        ret, frame = cap.read()
        if not ret:
            break
        frame_path = os.path.join(output_folder, f'frame{count}.png')
        cv2.imwrite(frame_path, frame)
        count += 1

    cap.release()
    cv2.destroyAllWindows()

def crop_to_central_data(image_path,number):
    """
    Crop the image to the central data area based on hardcoded coordinates.

    :param image_path: Path to the image file.
    :return: Cropped Image object.
    """
    crop_coordinates = (125, 120, 900, 890) 
    with Image.open(image_path) as img:
        # Crop the image to the predefined coordinates
        # Extract the directory and base name without extension
        directory = os.path.dirname(image_path)
        base_name = os.path.basename(image_path)
        name, extension = os.path.splitext(base_name)

        # Create the new file name and save the cropped image
        new_file_name = f"{'image_'}{number}{extension}"
        cropped_img = img.crop(crop_coordinates)
        cropped_img.save(os.path.join(directory, new_file_name))

def Fourrier_Analysis(frame):
    # FFT2 to identify horizontal wavenumber vector 
    # ChatGPT wrote this code and I tested/adapted from there 

    # remove enough colmuns or rows to make it square
    if frame.shape[0] > frame.shape[1]:
        frame = frame[:frame.shape[1],:]
    elif frame.shape[0] < frame.shape[1]:
        frame = frame[:,:frame.shape[0]]

    cloud_pattern = frame
    

    # Generate a sample image (replace this with your cloud probability distribution)
    image_size = cloud_pattern.shape[0]

    # Apply Fourier Transform
    fft_result = fft2(cloud_pattern)

    # Shift zero frequency components to the center
    fft_result_shifted = fftshift(fft_result)

    # Calculate amplitude
    amplitude = np.abs(fft_result_shifted)


    # Create 2D wavenumber array
    kx = np.fft.fftshift(np.fft.fftfreq(image_size, d=1/image_size)) 
    ky = np.fft.fftshift(np.fft.fftfreq(image_size, d=1/image_size))
    kx, ky = np.meshgrid(kx, ky)

    # Calculate total wavenumber array
    wavenumbers = np.sqrt(kx**2 + ky**2)


    # Calculate the polar coordinates
    radius = np.sqrt(kx**2 + ky**2)
    theta = np.arctan2(ky, kx)

    # Average the amplitude where the TOTAL wavenumber is between 4 and 20
    min_wavenumber = 2
    max_wavenumber = 20
    in_k_range = (radius>min_wavenumber) & (radius<max_wavenumber)


    # Create theta bins
    num_bins = 36  # 10 degrees each bin 
    theta_bins = np.linspace(np.min(theta), np.max(theta), num_bins + 1)

    # Calculate the mean amplitude in each theta bin
    mean_amplitude = np.zeros(num_bins)
    for i in range(num_bins):
        in_bin = (theta >= theta_bins[i]) & (theta < theta_bins[i + 1])
        mean_amplitude[i] = np.mean(amplitude * in_k_range * in_bin)


    # -------- chatGPT fits it, although really just argmax and max-mean suffices 
    # Define the bisinusoidal fit function
    def bisinusoidal_fit(theta, amplitude, phase, offset):
        return amplitude * np.sin(2*theta + phase) + offset
    # Define the sinusoidal fit function
    def sinusoidal_fit(theta, amplitude, phase, offset):
        return amplitude * np.sin(theta + phase) + offset

    # Initial guess for the fit parameters
    initial_guess = [1.0, 0.0, 0.0]

    # Fit the sinusoidal curve to the azimuthal variation
    popt, _ = curve_fit(bisinusoidal_fit, theta_bins[:-1], mean_amplitude, p0=initial_guess)
    # ------------


    # Plot the original image, amplitude, the annulus, and the fitted sinusoidal curve
    plt.figure(figsize=(16, 4))

    plt.subplot(141)
    plt.pcolormesh(cloud_pattern, cmap='viridis')
    plt.title('Original Image')

    plt.subplot(142)
    plt.pcolormesh(np.log1p(amplitude), cmap='viridis')
    plt.title('Log Amplitude Spectrum')

    plt.subplot(143)
    plt.pcolormesh(kx,ky, np.log1p(amplitude), cmap='viridis')
    ax = plt.gca()
    ax.set_xlim([-30, 30]) # type: ignore
    ax.set_ylim([-30, 30]) # type: ignore
    plt.contour(kx,ky,radius, levels=[min_wavenumber, max_wavenumber], colors='r', linewidths=2)
    plt.title('Annulus (Wavenumbers 4-20)')

    plt.subplot(144)
    plt.plot(theta_bins[:-1] *180/3.1415 + 360/num_bins/2, mean_amplitude)
    plt.title('amplitude vs. angle in 4-20 wavenumber band')
    plt.plot(theta_bins[:-1] *180/3.1415 + 360/num_bins/2, sinusoidal_fit(2*theta_bins[:-1], *popt), 'r-', label='Sinusoidal Fit')

    plt.plot() 

    plt.tight_layout()
    plt.show()

    argmx = np.argmax(mean_amplitude)
    maxangle = theta_bins[ argmx ]
    maxangledeg = theta_bins[ argmx ] *180/3.1415 + 360/num_bins/2

    print('angle is ', argmx, maxangle, maxangledeg, maxangledeg+180)
    print('strength is ',np.max(mean_amplitude)/np.mean(mean_amplitude) )

    # Print the fitted parameters
    #amplitude_fit, phase_fit, offset_fit = popt
    #print(f"Amplitude of the Fit: {amplitude_fit}")
    #print(f"Phase of the Fit: {phase_fit*1800./3.142}")
    #print(f"Offset of the Fit: {offset_fit}")

def FFT_Analysis(frames, plot_option='yes', plot_result ='yes'):

    data = frames
    data_shape = data.shape
    data_w = np.copy(data)
    for dim in range(1, 3):  # This changes from range(3) to range(1, 3)
        border_size = int(data_shape[dim] * 0.1)  # 10% of the dimension size
        if border_size < 1:
            border_size = 1  # Ensure at least one point gets the window applied
        
        full_window = np.hanning(2 * border_size)  # Full window for both sides
        window = np.ones(data_shape[dim])  # Create a window array full of ones
        window[:border_size] = full_window[:border_size]
        window[-border_size:] = full_window[border_size:]

        # Reshape the window to match the data dimension
        if dim == 1:
            window = window[np.newaxis, :, np.newaxis]
        elif dim == 2:
            window = window[np.newaxis, np.newaxis, :]

        # Multiply the data with the window
        data_w *= window

    # Perform 3D FFT
    fft_data = np.fft.fftn(data_w)
    fft_data = np.fft.fftshift(fft_data)

    # Calculate kl unitlessly for each pixel
    k = np.fft.fftshift(np.fft.fftfreq(data_shape[2], d=1/data_shape[2])) 
    l = np.fft.fftshift(np.fft.fftfreq(data_shape[1], d=1/data_shape[1]))
    k, l = np.meshgrid(k, l)

    k = k/(data_shape[2])
    l = l/(data_shape[1])

    # k_real = k/2000
    # l_real = l/2000
    # Calculate total wavenumber array
    kl = np.sqrt(k**2 + l**2)

    # Frequency components as unitless (normalized index positions)
    f = np.fft.fftfreq(data_shape[0])
    f = np.fft.fftshift(f)

    # Calculate unitless phase speed
    kl[kl == 0] = np.inf  # Avoid division by zero
    c = np.abs(f[:, np.newaxis, np.newaxis]) / kl  # c is unitless
    kl[kl == np.inf] = 0


    # Convert f to cycles per second (Hz) and kl to cycles per meter
    #f_real = f / 1800*data.shape[0]
    # f_real = f / 1800
    # kl_real = np.sqrt(k_real**2 + l_real**2)

    # kl_real[kl_real == 0] = np.inf 
    # c_real = np.abs(f_real[:, np.newaxis, np.newaxis]) / kl_real
    # kl_real[kl_real == np.inf] = 0

    kl = np.repeat(kl[np.newaxis, :, :], data_shape[0], axis=0)
    # kl_real = np.repeat(kl_real[np.newaxis, :, :], data_shape[0], axis=0)
    f = np.repeat(f[:, np.newaxis], data_shape[1], axis=1)
    f = np.repeat(f[:, :, np.newaxis], data_shape[2], axis=2)

    # f_real = np.repeat(f_real[:, np.newaxis], data_shape[1], axis=1)
    # f_real = np.repeat(f_real[:, :, np.newaxis], data_shape[2], axis=2)

    # Calculate the maximum value for consistent color scaling
    vmax = np.max(kl[0, :, :])

    #kl_cutoff = (1 / 5) * np.max(kl) 
    kl_cutoff_inf, kl_cutoff_sup = 1/250, 1/5
    c_cutoff_inf,c_cutoff_sup = 5,50 # pixels per frame

    if plot_option == 'yes':

        fig, axs = plt.subplots(1, 2, figsize=(12, 6))

        # Plot for the original data
        pc0 = axs[0].pcolormesh(kl[0, :, :], cmap='viridis', vmin=0, vmax=vmax)
        axs[0].set_title("KL amplitude spectrum, before filtering")
        axs[0].set_xlabel("KL")
        axs[0].grid(True)
        fig.colorbar(pc0, ax=axs[0])  # Attach the colorbar to the first plot

        # Apply filtering
        kl_disp= kl.copy()
        #filter_mask = (kl <=  kl_cutoff ) 
        filter_mask = (kl >= kl_cutoff_inf) & (kl <= kl_cutoff_sup )
        kl_disp[~filter_mask] = 0

        # Plot for the filtered data
        pc1 = axs[1].pcolormesh(kl_disp[0, :, :], cmap='viridis', vmin=0, vmax=vmax)
        axs[1].set_title("KL amplitude spectrum, after filtering")
        axs[1].set_xlabel("KL")
        axs[1].grid(True)
        fig.colorbar(pc1, ax=axs[1])  # Attach the colorbar to the second plot

        plt.show()

        ## Now the same for Phase Speed
        # Calculate the maximum value for consistent color scaling
        vmax = np.percentile(c[0, :, :],99)

        fig, axs = plt.subplots(1, 2, figsize=(12, 6))

        # Plot for the original data
        pc0 = axs[0].pcolormesh(c[0, :, :], cmap='viridis', vmin=0, vmax=vmax)
        axs[0].set_title("Phase speed amplitude spectrum, before filtering")
        axs[0].set_xlabel("Phase speed")
        axs[0].grid(True)
        fig.colorbar(pc0, ax=axs[0])  # Attach the colorbar to the first plot

        # Apply filtering
        c_disp= c.copy()
        filter_mask = (c > c_cutoff_inf) & (c < c_cutoff_sup)
        c_disp[~filter_mask] = 0

        # Plot for the filtered data
        pc1 = axs[1].pcolormesh(c_disp[0, :, :], cmap='viridis', vmin=0, vmax=vmax)
        axs[1].set_title("Phase speed amplitude spectrum, after filtering")
        axs[1].set_xlabel("Phase speed")
        axs[1].grid(True)
        fig.colorbar(pc1, ax=axs[1])  # Attach the colorbar to the second plot

        plt.show()

    # Define and apply a filter based on unitless criteria
    # filter_mask = (kl_real < 1/10000) & (kl_real > 1/100000) & (np.abs(f) > (1 / 24.0)) & (c_real > 5) & (c_real < 50) 
    filter_mask = (c > c_cutoff_inf) & (c < c_cutoff_sup) & (kl >= kl_cutoff_inf) & (kl <= kl_cutoff_sup )

    fft_data[~filter_mask] = 0

    if plot_option == 'yes':

        # Plotting the distribution of kl and c values 
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))
        axs[0].hist(kl.flatten(), bins=500, color='blue', alpha=0.7)
        axs[0].set_title("wavenumber cycles per pixel")
        axs[0].set_xlabel("wavenumber (cpm)")
        axs[0].set_ylabel("Frequency")
        axs[0].grid(True)
        filter_mask = (kl >= kl_cutoff_inf) & (kl <= kl_cutoff_sup )
        axs[0].axvline(x=kl_cutoff_inf, color='black', linestyle='--')
        axs[0].axvline(x=kl_cutoff_sup, color='black', linestyle='--')
        # overlay an histogram of the kl values with the filter mask applied
        axs[0].hist(kl[filter_mask].flatten(), bins=500, color='green', alpha=0.7, range = (0, np.max(kl)))
        axs[0].grid(True)

        axs[1].hist(c.flatten(), bins=500, color='red', alpha=0.7, range = (0, np.percentile(c,99.99)), log=True)
        axs[1].set_title("Phase Speed (Pixel per Frame)")
        axs[1].set_xlabel("Phase Speed (Pixel per Frame)")
        axs[1].set_ylabel("Frequency")
        axs[1].axvline(x=c_cutoff_inf, color='black', linestyle='--')
        axs[1].axvline(x=c_cutoff_sup, color='black', linestyle='--')
        filter_mask = (c > c_cutoff_inf) & (c < c_cutoff_sup) 
        # overlay an histogram of the c values with the filter mask applied
        axs[1].hist(c[filter_mask].flatten(), bins=500, color='green', alpha=0.7, range = (0, np.percentile(c,99.99)), log=True)
        axs[1].grid(True)


    # Inverse FFT and return
    fft_data = np.fft.ifftshift(fft_data)
    filtered_data = np.fft.ifftn(fft_data)
    filtered_data = np.real(filtered_data)

    if plot_result == 'yes':
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))
        axs[0].imshow(data[0], cmap='gray')
        axs[0].set_title("First Frame")
        #fig.colorbar(axs[0].imshow(data[0], cmap='gray'), ax=axs[0])
        axs[1].imshow(filtered_data[0], cmap='gray')
        axs[1].set_title("First Filtered Frame")
        #fig.colorbar(axs[1].imshow(filtered_data[0], cmap='gray'), ax=axs[1])
        plt.show()

    return filtered_data

def PCRegression(data1, data2, plot_option='no', ax=None):
    # Example datasets
    x = data1.flatten()
    y = data2.flatten()

    # Prepare and standardize data
    data = np.column_stack((x, y))
    data -= np.mean(data, axis=0)

    # Apply PCA
    pca = PCA(n_components=2)
    pca.fit(data)

    # Get the principal components
    components = pca.transform(data)

    # Use the first principal component for regression
    pc1 = components[:, 0].reshape(-1, 1)

    # Linear regression
    model = LinearRegression()
    model.fit(pc1, y)

    if plot_option == 'yes':
        if ax is None:
            ax = plt.gca()  # Get current axis if none is provided

        # Plotting using the provided axis
        #ax.scatter(pc1, y, alpha=0.5, label='Data Points')
        #ax.plot(pc1, model.predict(pc1), color='red', label='Regression Line')
        ax.scatter(x, y, alpha=0.5, label='Data Points')
        ax.plot(-1*pc1, model.predict(pc1), color='red', label='Regression Line')
        ax.set_title('Principal Component Regression')
        ax.legend()

        # Optionally, you can still show the plot if the function is responsible for it
        if ax is plt.gca():  # Checks if the default axis was used
            plt.show()

        # Print variance explained, if needed
        print(f"Variance explained by the first component: {pca.explained_variance_ratio_[0]}")

def Get_Candidates(x,y,threshold,plot_option='yes'):
    # x is brightness, y is conv

    # Create a structured array to store the candidates
    candidates = np.zeros((y.shape[0], y.shape[1]), dtype=[('con', float), ('bri', float)])
    candidates['con'] = y
    candidates['bri'] = x

    # Calculate the line equation y = mx where m = std(y) / std(x)
    m = np.std(y) / np.std(x)
    b = np.mean(y) - m * np.mean(x)

    # Calculate perpendicular distances to the line y = mx + b
    distances = np.abs(y - (m * x + b)) / np.sqrt(m**2 + 1)

    # Set threshold for closeness (e.g., within 50 units)
    close_to_line = distances < threshold

    if plot_option == 'yes':
        plt.figure(figsize=(10, 5))
        # plt.scatter(x, y, alpha=0.5, label='All Points')
        # plt.scatter(*close_points, color='red', label='Close to 45-degree line')
        plt.scatter(candidates['bri'][close_to_line], candidates['con'][close_to_line], s=1, color='red')
        plt.scatter(candidates['bri'][~close_to_line], candidates['con'][~close_to_line], s=1, color='blue')
        plt.plot(np.array([x.min(), x.max()]), np.array([m*x.min() + b, m*x.max() + b]), 'g--', label='Adjusted 45-degree line')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Points Close to the Adjusted 45-Degree Line')
        plt.legend(loc='upper left')
        plt.show()

    # Retrieve the coordinates of the points close to the 1:1 line
    close_points_indices = np.nonzero(close_to_line)
    close_points_coordinates = list(zip(close_points_indices[0], close_points_indices[1]))

    if plot_option == 'yes':
        # calculate percentage of points close to the line
        total_points = candidates.shape[0] * candidates.shape[1]
        percentage = 100 * len(close_points_coordinates) / total_points
        print(f'Percentage of points close to the line: {percentage:.2f}%')

    return close_points_coordinates

def Plot_Candidates (frames_data, frames_index, close_points_coordinates):
    
    plt.figure(figsize=(10, 5))
    plt.imshow(frames_data[frames_index[0]], cmap='gray')
    plt.scatter([p[1] for p in close_points_coordinates], [p[0] for p in close_points_coordinates], s=1, color='red')
    plt.title('Potential candidates for GW')
    plt.show()

def Show_Power_Spectrum(frames,index=0,theta_coeff=0,spat_angle_cutoff_inf= -np.pi/2, spat_angle_cutoff_sup= np.pi/2, prop_angle_cutoff_inf= -np.pi/2, prop_angle_cutoff_sup= np.pi/2,angle_diff=0.1,
                        kl_cutoff_inf= 1/250, kl_cutoff_sup= 1/5 , c_cutoff_inf= 5,c_cutoff_sup= 50):
    data = frames
    data_shape = data.shape
    data_w = np.copy(data)
    for dim in range(1, 3):  # This changes from range(3) to range(1, 3)
        border_size = int(data_shape[dim] * 0.1)  # 10% of the dimension size
        if border_size < 1:
            border_size = 1  # Ensure at least one point gets the window applied
        
        full_window = np.hanning(2 * border_size)  # Full window for both sides
        window = np.ones(data_shape[dim])  # Create a window array full of ones
        window[:border_size] = full_window[:border_size]
        window[-border_size:] = full_window[border_size:]

        # Reshape the window to match the data dimension
        if dim == 1:
            window = window[np.newaxis, :, np.newaxis]
        elif dim == 2:
            window = window[np.newaxis, np.newaxis, :]

        # Multiply the data with the window
        data_w *= window

    # Perform 3D FFT
    fft_data = np.fft.fftn(data_w)
    fft_data = np.fft.fftshift(fft_data)

    # Calculate kl unitlessly for each pixel
    k = np.fft.fftshift(np.fft.fftfreq(data_shape[2], d=1/data_shape[2])) 
    l = np.fft.fftshift(np.fft.fftfreq(data_shape[1], d=1/data_shape[1]))
    k, l = np.meshgrid(k, l)

    k = k/(data_shape[2])
    l = l/(data_shape[1])

    spat_angle = np.arctan2(-l, k) 
    spat_angle = np.repeat(spat_angle[np.newaxis, :, :], data_shape[0], axis=0)
    prop_angle = np.arctan2(-k, l)
    prop_angle = np.repeat(prop_angle[np.newaxis, :, :], data_shape[0], axis=0)
    
    angle_difference = np.abs(prop_angle - spat_angle- np.pi/2) -np.pi

    # Calculate total wavenumber array
    kl = np.sqrt(k**2 + l**2)

    # Frequency components as unitless (normalized index positions)
    f = np.fft.fftfreq(data_shape[0])
    f = np.fft.fftshift(f)

    # Calculate unitless phase speed
    kl[kl == 0] = np.inf  # Avoid division by zero
    c = np.abs(f[:, np.newaxis, np.newaxis]) / kl  # c is unitless
    kl[kl == np.inf] = 0

    kl = np.repeat(kl[np.newaxis, :, :], data_shape[0], axis=0)
    f = np.repeat(f[:, np.newaxis], data_shape[1], axis=1)
    f = np.repeat(f[:, :, np.newaxis], data_shape[2], axis=2)

    # Calculate the polar coordinates
    radius = np.sqrt(k**2 + l**2)
    theta = np.arctan2(k, l)

    # Define and apply a filter based on unitless criteria
    # filter_mask = (np.abs(f) > (1 / 24.0)) 
    filter_mask = (c > c_cutoff_inf) & (c < c_cutoff_sup) & (spat_angle >= spat_angle_cutoff_inf) & (spat_angle <= spat_angle_cutoff_sup ) & (prop_angle >= prop_angle_cutoff_inf ) & (prop_angle <= prop_angle_cutoff_sup) & (kl >= kl_cutoff_inf) & (kl <= kl_cutoff_sup ) & (angle_difference <= angle_diff) & (angle_difference >= -angle_diff)

    fft_data = fft_data*(1+theta_coeff*np.cos(2*theta))
    data_power = np.abs(fft_data)
    fft_data[~filter_mask] = 0
    data_power_filt = np.abs(fft_data)

    log_data_power = np.log1p(data_power)
    log_data_power_filt = np.log1p(data_power_filt)

    # Set a stringent threshold based on the power distribution
    threshold_value = np.percentile(log_data_power[index,:,:], 80)  # Adjust percentile as needed
    mask = log_data_power[index,:,:] > threshold_value

    # Find contours on the thresholded image
    contours = measure.find_contours(mask, 0.8)  # Adjust level as needed based on visualization

    # Compute the centroid of the high-power region
    mask_indices = np.where(mask)
    power_centroid = [np.mean(mask_indices[0]), np.mean(mask_indices[1])]

    # Find the contour that is closest to the centroid of the high-power region
    min_distance = float('inf')
    best_contour = None
    for contour in contours:
        contour_centroid = [np.mean(contour[:, 0]), np.mean(contour[:, 1])]
        distance = np.linalg.norm(np.array(contour_centroid) - np.array(power_centroid))
        if distance < min_distance:
            min_distance = distance
            best_contour = contour

    # Extract x and y coordinates
    X = best_contour[:, 1:2]  # x coordinates # type: ignore
    Y = best_contour[:, 0:1]  # y coordinates # type: ignore

    # Formulate and solve the least squares problem ||Ax - b ||^2
    A = np.hstack([X**2, X * Y, Y**2, X, Y])
    b = np.ones_like(X)
    x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)  # Use rcond=None for future compatibility
    x = x.squeeze()

    # Print the equation of the ellipse in standard form
    print('The ellipse is given by {:.3f}x^2 + {:.3f}xy + {:.3f}y^2 + {:.3f}x + {:.3f}y = 1'.format(x[0], x[1], x[2], x[3], x[4]))

    # Plot the fitted ellipse
    x_coord = np.linspace(np.min(X), np.max(X), 300)
    y_coord = np.linspace(np.min(Y), np.max(Y), 300)
    X_coord, Y_coord = np.meshgrid(x_coord, y_coord)
    Z_coord = x[0] * X_coord**2 + x[1] * X_coord * Y_coord + x[2] * Y_coord**2 + x[3] * X_coord + x[4] * Y_coord
    fig, ax = plt.subplots()
    ax.pcolormesh(log_data_power[index,:,:], shading='auto')
    ax.plot(best_contour[:, 1], best_contour[:, 0], linewidth=2, color='blue')  # type: ignore # Plot the best contour
    ax.contour(X_coord, Y_coord, Z_coord, levels=[1], colors='red', linewidths=2)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.show()


    plt.show()
    if theta_coeff != 0:
        # Calculate amplitude
        amplitude = np.abs(fft_data[index,:,:])
        amplitude_tilted = np.abs(fft_data[index,:,:]*(1+theta_coeff*np.cos(2*theta)))

        frame_tilted = np.real(np.fft.ifft2(np.fft.ifftshift(fft_data)))

        plt.figure(figsize=(16*2, 4*2))

        plt.subplot(141)
        plt.pcolormesh(k,l,np.log1p(amplitude), cmap='viridis')
        plt.title('Log Amplitude Spectrum')
        plt.colorbar()

        plt.subplot(142)
        plt.pcolormesh(k,l,np.log1p(amplitude_tilted), cmap='viridis')
        plt.title('Tilted Log Amplitude Spectrum')
        plt.colorbar()

        plt.tight_layout()
        plt.show()

        plt.figure(figsize=(16*2, 4*2))

        plt.subplot(141)
        plt.imshow(frames[index,:,:], cmap='gray')
        plt.title('Original Image')

        plt.subplot(142)
        plt.imshow(frame_tilted[index,:,:], cmap='gray')
        plt.title('Tilted Image')

def bisinusoidal_func(theta, A1, phi1, C):
    return A1 * np.cos(2 * theta + phi1) +  C

def analyze_frame_SFilt(frames, frame_index, square_size=256, overlap=0.5, strength=0.75, kl_cutoff_inf=0.01, 
                  kl_cutoff_sup=0.1, wavenumber_step=0.01, c_cutoff_inf=0, c_cutoff_sup=500, Filter=True, log=True, plot_condition=False):
    
    data_shape = frames.shape
    all_directions = []
    all_amplitudes = []


    step_size = int(square_size * (1 - overlap))
    wavenumber_values = np.arange(kl_cutoff_inf, kl_cutoff_sup, wavenumber_step)

    k = fftshift(np.fft.fftfreq(data_shape[1], d=1/data_shape[1]))
    l = fftshift(np.fft.fftfreq(data_shape[0], d=1/data_shape[0]))
    k, l = np.meshgrid(k, l)
    k = k / data_shape[1]
    l = l / data_shape[0]

    radius = np.sqrt(k**2 + l**2)
    theta = np.arctan2(l, k)

    kl = np.repeat(radius[:, :,np.newaxis], data_shape[2], axis=2)

    # Frequency components as unitless (normalized index positions)
    f = np.fft.fftfreq(data_shape[2])
    f = np.repeat(f[np.newaxis, :], data_shape[1], axis=0)
    f = np.repeat(f[np.newaxis, :, :], data_shape[0], axis=0)
    f = np.fft.fftshift(f)

    # Calculate unitless phase speed
    kl[kl == 0] = np.inf  # Avoid division by zero
    c = np.abs(f) / kl  # c is unitless
    kl[kl == np.inf] = 0


    if Filter:
        filter_mask = (c >= c_cutoff_inf) & (c <= c_cutoff_sup) & (kl >= kl_cutoff_inf) & (kl <= kl_cutoff_sup) 

        # Perform FFT on the entire frame
        fft_data = np.fft.fftn(frames)
        fft_data = np.fft.fftshift(fft_data)

        # Apply the mask to the FFT data
        fft_data[~filter_mask] = 0

        # Perform IFFT on the filtered data
        fft_data = np.fft.ifftshift(fft_data)
        fft_data = np.fft.ifftn(fft_data)
        frames = np.real(fft_data)
        
    frame = frames[:,:,frame_index]
    data_shape = frame.shape

    for i in range(0, data_shape[0] - square_size + 1, step_size):
        for j in range(0, data_shape[1] - square_size + 1, step_size):
            # Create a mask for the entire frame
            mask = np.zeros(data_shape)

            # Create a 2D radial gradient mask for the area outside the square
            y, x = np.ogrid[:data_shape[0], :data_shape[1]]
            center_y, center_x = i + square_size // 2, j + square_size // 2
            distance_from_center = np.sqrt((x - center_x)**2 + (y - center_y)**2)
            max_distance = strength * np.sqrt((data_shape[0] / 2)**2 + (data_shape[1] / 2)**2)
            gradient_mask = 1 - np.clip(distance_from_center / max_distance, 0, 1)

            # Apply the gradient mask to the area outside the central square
            mask[i:i+square_size, j:j+square_size] = 1
            mask = np.maximum(mask, gradient_mask)

            # Apply the mask to the frame
            blurred_frame = frame * mask

            # Perform FFT on the entire frame
            fft_data = fftn(blurred_frame)
            fft_data = fftshift(fft_data) 

            amplitude = np.abs(fft_data)

            directions = []
            amplitudes = []

            for wavenumber in wavenumber_values:
                mask = (radius > (wavenumber - wavenumber_step/2)) & (radius < (wavenumber + wavenumber_step/2))
                if np.any(mask):
                    extracted_data = amplitude[mask]
                    theta_flat = theta[mask]
                    
                    # Create theta bins
                    num_bins = 36  # 10 degrees each bin
                    theta_bins = np.linspace(np.min(theta_flat), np.max(theta_flat), num_bins + 1)
                    mean_amplitude = np.zeros(num_bins)

                    for k in range(num_bins):
                        in_bin = (theta_flat >= theta_bins[k]) & (theta_flat < (theta_bins[k + 1] if k < num_bins - 1 else theta_bins[k]))
                        bin_data = extracted_data[in_bin]
                        if len(bin_data) > 0:
                            mean_amplitude[k] = np.mean(bin_data)
                        else:
                            mean_amplitude[k] = 0  # or another default value

                    # Perform bisinusoidal fitting
                    initial_guess = [np.max(mean_amplitude), 0, np.mean(mean_amplitude)]
                    popt, _ = curve_fit(bisinusoidal_func, theta_bins[:-1], mean_amplitude, p0=initial_guess)

                    # Calculate direction and amplitude based on the fit
                    theta_fit = theta_bins[:-1]
                    fitted_amplitude = bisinusoidal_func(theta_fit, *popt)
                    max_indices = np.argsort(fitted_amplitude)[-2:]  # Indices of the two largest values

                    direction_1 = theta_fit[max_indices[0]] * 180 / np.pi   # X value for the first maximum Y in degrees
                    amplitude_1 = np.max([fitted_amplitude[max_indices[0]],fitted_amplitude[max_indices[1]]])  # First maximum Y value

                    direction_2 = theta_fit[max_indices[1]] * 180 / np.pi   # X value for the second maximum Y in degrees

                    # Plot the amplitude vs. angle and the fit if Plot is True
                    if plot_condition:
                        plt.figure(figsize=(10, 6))
                        plt.plot(theta_bins[:-1] * 180 / np.pi, mean_amplitude, label='Mean Amplitude')
                        plt.plot(theta_bins[:-1] * 180 / np.pi, bisinusoidal_func(theta_bins[:-1], *popt), 'r-', label='Bisinusoidal Fit')
                        plt.title(f'Amplitude vs. Angle for Wavenumber {wavenumber}')
                        plt.xlabel('Angle (degrees)')
                        plt.ylabel('Amplitude')
                        plt.legend()
                        plt.grid(True)
                        plt.show()
                        print(f'Fitted parameters for wavenumber {wavenumber}: A1 = {popt[0]}, phi1 = {popt[1]}, C = {popt[2]}')

                    # Print the calculated directions and amplitudes
                    # print(f'Calculated directions: {direction_1} degrees, {direction_2} degrees')
                    # print(f'Calculated amplitude: {amplitude_1}')

                    # Store direction and amplitude data
                    directions.append((direction_1, direction_2))
                    amplitudes.append(amplitude_1)

            all_directions.append(directions)
            all_amplitudes.append(amplitudes)

    return all_directions, all_amplitudes, wavenumber_values

def visualize_results_SFilt(frames, frame_index, all_directions, all_amplitudes, wavenumber_values, 
                      square_number, square_size=256, overlap=0.5, strength=0.75,
                       kl_cutoff_inf=0.01, kl_cutoff_sup=0.1, c_cutoff_inf=0, 
                       c_cutoff_sup=500, Filter=True,Amp=False, log=True, Km=True):
    
    data_shape = frames.shape

    # Filter the frame first
    if Filter:
        # Perform FFT on the entire frame
        fft_data = np.fft.fftn(frames)
        fft_data = np.fft.fftshift(fft_data)

        # Calculate kl unitlessly for each pixel
        k = np.fft.fftshift(np.fft.fftfreq(data_shape[1], d=1/data_shape[1])) 
        l = np.fft.fftshift(np.fft.fftfreq(data_shape[0], d=1/data_shape[0]))
        k, l = np.meshgrid(k, l)

        k = k/(data_shape[1])
        l = l/(data_shape[0])


        kl = np.sqrt(k**2 + l**2)
        kl = np.repeat(kl[:, :,np.newaxis], data_shape[2], axis=2)


        # Frequency components as unitless (normalized index positions)
        f = np.fft.fftfreq(data_shape[2])
        f = np.repeat(f[np.newaxis, :], data_shape[1], axis=0)
        f = np.repeat(f[np.newaxis, :, :], data_shape[0], axis=0)
        f = np.fft.fftshift(f)

        # Calculate unitless phase speed
        kl[kl == 0] = np.inf  # Avoid division by zero
        c = np.abs(f) / kl  # c is unitless
        kl[kl == np.inf] = 0


        filter_mask = (c >= c_cutoff_inf) & (c <= c_cutoff_sup) & (kl >= kl_cutoff_inf) & (kl <= kl_cutoff_sup) 

        # Apply the mask to the FFT data

        fft_data[~filter_mask] = 0

        fft_data = np.fft.ifftshift(fft_data)
        fft_data = np.fft.ifftn(fft_data)
        frames = np.real(fft_data)

        if  Amp == True:
            frames = np.abs(frames)
        
    frame = frames[:,:,frame_index]
    data_shape = frame.shape

    step_size = int(square_size * (1 - overlap))
    num_squares_per_row = (data_shape[1] - square_size) // step_size + 1
    
    i = (square_number // num_squares_per_row) * step_size
    j = (square_number % num_squares_per_row) * step_size

    # Create a mask for the entire frame
    mask = np.zeros(data_shape)

    # Create a 2D radial gradient mask for the area outside the square
    y, x = np.ogrid[:data_shape[0], :data_shape[1]]
    center_y, center_x = i + square_size // 2, j + square_size // 2
    distance_from_center = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    max_distance = strength * np.sqrt((data_shape[0] / 2)**2 + (data_shape[1] / 2)**2)
    gradient_mask = 1 - np.clip(distance_from_center / max_distance, 0, 1)

    # Apply the gradient mask to the area outside the central square
    mask[i:i+square_size, j:j+square_size] = 1
    mask = np.maximum(mask, gradient_mask)

    # Apply the mask to the frame
    blurred_frame = frame * mask

    # Plot the blurred frame and the corresponding rose plot
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    ax_blur = axs[0]
    ax_polar = fig.add_subplot(122, polar=True)  # Ensure the second subplot is a polar plot
    cmap = plt.get_cmap('coolwarm')
    colors = cmap(np.linspace(0, 1, len(wavenumber_values)))

    directions = all_directions[square_number]
    amplitudes = all_amplitudes[square_number]

    for i, wavenumber in enumerate(wavenumber_values):
        if i < len(directions) and i < len(amplitudes):
            direction = directions[i]
            amplitude = amplitudes[i]

            if log:
                # Apply log scale to amplitudes
                amplitude = np.log(amplitude + 1)  # Adding 1 to avoid log(0)

            # Wrap angles to the range [0, 360] degrees
            direction_wrapped = np.mod(direction, 360)

            # Plot the histogram for both directions with the same amplitude
            for d in direction_wrapped:
                if Km == True:
                    ax_polar.bar(np.deg2rad(d), amplitude, width=np.deg2rad(10), color=colors[i], alpha=0.6, edgecolor='k', label=f'{2*1/wavenumber:.3f}' if d == direction_wrapped[0] else "")
                else:
                    ax_polar.bar(np.deg2rad(d), amplitude, width=np.deg2rad(10), color=colors[i], alpha=0.6, edgecolor='k', label=f'{2*1/wavenumber:.3f}' if d == direction_wrapped[0] else "")
                
    ax_polar.set_title('Circular Spectrum of Directional Angles')
    ax_polar.set_theta_zero_location('N')
    ax_polar.set_theta_direction(-1)
    ax_polar.legend(loc='upper right', bbox_to_anchor=(1.1, 1.1))

    ax_blur.pcolormesh(blurred_frame, cmap='gray')
    ax_blur.set_title(f'Square {square_number} Blurred Frame')

    plt.show()
   
def Visualize_Filter(frames,frame,kl_cutoff_inf=0.01, kl_cutoff_sup=0.1, c_cutoff_inf=0, c_cutoff_sup=500, Amp=False):
    data = frames
    data_shape = data.shape
    data_w = np.copy(data)

    # Needs a rework for dynamic determination of the time dimension
    # for dim in range(1, 3):  # This changes from range(3) to range(1, 3)
    #     border_size = int(data_shape[dim] * 0.1)  # 10% of the dimension size
    #     if border_size < 1:
    #         border_size = 1  # Ensure at least one point gets the window applied
        
    #     full_window = np.hanning(2 * border_size)  # Full window for both sides
    #     window = np.ones(data_shape[dim])  # Create a window array full of ones
    #     window[:border_size] = full_window[:border_size]
    #     window[-border_size:] = full_window[border_size:]

    #     # Reshape the window to match the data dimension
    #     if dim == 1:
    #         window = window[np.newaxis, :, np.newaxis]
    #     elif dim == 2:
    #         window = window[np.newaxis, np.newaxis, :]

    #     # Multiply the data with the window
    #     data_w *= window

    # Perform 3D FFT
    # fft_data = np.fft.fftn(data_w)
    fft_data = np.fft.fftn(data)
    fft_data = np.fft.fftshift(fft_data)

    # Calculate kl unitlessly for each pixel
    k = np.fft.fftshift(np.fft.fftfreq(data_shape[1], d=1/data_shape[1])) 
    l = np.fft.fftshift(np.fft.fftfreq(data_shape[0], d=1/data_shape[0]))
    k, l = np.meshgrid(k, l)

    k = k/(data_shape[1])
    l = l/(data_shape[0])


    kl = np.sqrt(k**2 + l**2)
    kl = np.repeat(kl[:, :,np.newaxis], data_shape[2], axis=2)


    # Frequency components as unitless (normalized index positions)
    f = np.fft.fftfreq(data_shape[2])
    f = np.repeat(f[np.newaxis, :], data_shape[1], axis=0)
    f = np.repeat(f[np.newaxis, :, :], data_shape[0], axis=0)
    f = np.fft.fftshift(f)

    # Calculate unitless phase speed
    kl[kl == 0] = np.inf  # Avoid division by zero
    c = np.abs(f) / kl  # c is unitless
    kl[kl == np.inf] = 0


    filter_mask = (c >= c_cutoff_inf) & (c <= c_cutoff_sup) & (kl >= kl_cutoff_inf) & (kl <= kl_cutoff_sup) 

    # Apply the mask to the FFT data

    fft_data[~filter_mask] = 0

    fft_data = np.fft.ifftshift(fft_data)
    fft_data = np.fft.ifftn(fft_data)
    frame = np.real(fft_data)

    if  Amp == True:
        frame = np.abs(frame)

    data_shape = frame.shape

    plt.figure(figsize=(16*2, 4*2))

    plt.subplot(141)
    plt.pcolormesh(frames[:,:,0], cmap='gray')
    plt.title('Original Image')

    plt.subplot(142)
    plt.pcolormesh(frame[:,:,0], cmap='gray')
    plt.title('Filtered Image')

    plt.tight_layout()
    plt.show()   

def Visualize_Filter_3D(frames,kl_inf = 0.00666666666, kl_sup = 0.06666666666, c_inf = 12, c_sup = 17):
        
    data_shape = frames.shape

    k = fftshift(np.fft.fftfreq(data_shape[1], d=1/data_shape[1]))
    l = fftshift(np.fft.fftfreq(data_shape[0], d=1/data_shape[0]))
    k, l = np.meshgrid(k, l)
    k = k / data_shape[1]
    l = l / data_shape[0]

    radius = np.sqrt(k**2 + l**2)
    theta = np.arctan2(l, k)

    kl = np.repeat(radius[:, :,np.newaxis], data_shape[2], axis=2)

    # Frequency components as unitless (normalized index positions)
    f = np.fft.fftfreq(data_shape[2])
    f = np.repeat(f[np.newaxis, :], data_shape[1], axis=0)
    f = np.repeat(f[np.newaxis, :, :], data_shape[0], axis=0)
    f = np.fft.fftshift(f)

    # Calculate unitless phase speed
    kl[kl == 0] = np.inf  # Avoid division by zero
    c = np.abs(f) / kl  # c is unitless
    kl[kl == np.inf] = 0

    filter_mask = (c >= c_inf) & (c <= c_sup) & (kl >= kl_inf) & (kl <= kl_sup)

    # Create a 3D plot
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Get the coordinates of the True values
    x, y, z = np.nonzero(filter_mask)

    # Plot the True values
    ax.scatter(x, y, z, c='b', marker='o', alpha=0.1)


    # Set labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Set limits to better visualize
    ax.set_xlim(filter_mask.shape[0]/4, 3*filter_mask.shape[0]/4)
    ax.set_ylim(filter_mask.shape[1]/4, 3*filter_mask.shape[1]/4)
    ax.set_zlim(0, filter_mask.shape[2])

    # Adjust the viewing angle
    ax.view_init(elev=20, azim=45)  # Adjust elev and azim to change the angle

    plt.show()

def compute_phase_shift(frame1, frame2):
    epsilon = 1e-10  # Small value to avoid division by zero
    fft_frame1 = fft2(frame1)
    fft_frame2 = fft2(frame2)
    cross_power_spectrum = (fft_frame1 * np.conj(fft_frame2)) / (np.abs(fft_frame1 * np.conj(fft_frame2)) + epsilon)
    phase_corr = np.abs(ifft2(cross_power_spectrum))
    max_corr_idx = np.unravel_index(np.argmax(phase_corr), phase_corr.shape)
    shift_vector = np.array(max_corr_idx) - np.array(frame1.shape) // 2
    return shift_vector

def filter_frames(frames, kl_cutoff_inf=0.01, kl_cutoff_sup=0.1, c_cutoff_inf=0, c_cutoff_sup=500):
    data_shape = frames.shape
    fft_data = np.fft.fftn(frames)
    fft_data = np.fft.fftshift(fft_data)

    k = np.fft.fftshift(np.fft.fftfreq(data_shape[1], d=1/data_shape[1]))
    l = np.fft.fftshift(np.fft.fftfreq(data_shape[0], d=1/data_shape[0]))
    k, l = np.meshgrid(k, l)
    k = k / data_shape[1]
    l = l / data_shape[0]

    radius = np.sqrt(k**2 + l**2)

    kl = np.repeat(radius[:, :, np.newaxis], data_shape[2], axis=2)

    f = np.fft.fftfreq(data_shape[2])
    f = np.repeat(f[np.newaxis, :], data_shape[1], axis=0)
    f = np.repeat(f[np.newaxis, :, :], data_shape[0], axis=0)
    f = np.fft.fftshift(f)

    kl[kl == 0] = np.inf
    c = np.abs(f) / kl
    kl[kl == np.inf] = 0

    filter_mask = (c >= c_cutoff_inf) & (c <= c_cutoff_sup) & (kl >= kl_cutoff_inf) & (kl <= kl_cutoff_sup)
    fft_data[~filter_mask] = 0

    fft_data = np.fft.ifftshift(fft_data)
    filtered_frames = np.real(np.fft.ifftn(fft_data))

    return filtered_frames

def infer_propagation_direction(frames, square_number=None, square_size=256, strength=0.5, overlap=0.5, kl_cutoff_inf=0.01, kl_cutoff_sup=0.1, c_cutoff_inf=0, c_cutoff_sup=500, Filter=True):
    if Filter:
        frames = filter_frames(frames, kl_cutoff_inf, kl_cutoff_sup, c_cutoff_inf, c_cutoff_sup)
    
    data_shape = frames.shape
    mask = None
    if square_number is not None:
        step_size = int(square_size * (1 - overlap))
        num_squares_per_row = (data_shape[1] - square_size) // step_size + 1
        
        i = (square_number // num_squares_per_row) * step_size
        j = (square_number % num_squares_per_row) * step_size

        # Create a mask for the entire frame
        mask = np.zeros((data_shape[0], data_shape[1]))

        # Create a 2D radial gradient mask for the area outside the square
        y, x = np.ogrid[:data_shape[0], :data_shape[1]]
        center_y, center_x = i + square_size // 2, j + square_size // 2
        distance_from_center = np.sqrt((x - center_x)**2 + (y - center_y)**2)
        max_distance = strength * np.sqrt((data_shape[0] / 2)**2 + (data_shape[1] / 2)**2)
        gradient_mask = 1 - np.clip(distance_from_center / max_distance, 0, 1)

        # Apply the gradient mask to the area outside the central square
        mask[i:i+square_size, j:j+square_size] = 1
        mask = np.maximum(mask, gradient_mask)

        # Apply the mask to all frames
        for frame_index in range(data_shape[2]):
            frames[:, :, frame_index] *= mask
    
    num_frames = frames.shape[2]
    total_shift = np.array([0, 0], dtype=float)
    for frame_index in range(num_frames - 1):
        frame1 = frames[:, :, frame_index]
        frame2 = frames[:, :, frame_index + 1]
        shift_vector = compute_phase_shift(frame1, frame2)
        total_shift += shift_vector
    avg_shift = total_shift / (num_frames - 1)
    avg_direction = np.arctan2(avg_shift[0], avg_shift[1])
    # Adjust the direction to make 0 degrees at the north and measure clockwise
    avg_direction = (np.pi / 2) - avg_direction
    # Ensure the angle is within the range [0, 360) degrees
    avg_direction_deg = (np.degrees(avg_direction) + 360) % 360
    return avg_direction_deg, avg_shift, mask, (i, j, square_size)

def analyze_amplitude_directions(blurred_frame, radius, theta, wavenumber, wavenumber_step, num_bins=36):
    mask = (radius > (wavenumber - wavenumber_step/2)) & (radius < (wavenumber + wavenumber_step/2))
    if np.any(mask):
        extracted_data = blurred_frame[mask]
        theta_flat = theta[mask]
        
        theta_bins = np.linspace(np.min(theta_flat), np.max(theta_flat), num_bins + 1)
        mean_amplitude = np.zeros(num_bins)

        for k in range(num_bins):
            in_bin = (theta_flat >= theta_bins[k]) & (theta_flat < (theta_bins[k + 1] if k < num_bins - 1 else theta_bins[k]))
            bin_data = extracted_data[in_bin]
            if len(bin_data) > 0:
                mean_amplitude[k] = np.mean(bin_data)
            else:
                mean_amplitude[k] = 0

        initial_guess = [np.max(mean_amplitude), 0, np.mean(mean_amplitude)]
        popt, _ = curve_fit(bisinusoidal_func, theta_bins[:-1], mean_amplitude, p0=initial_guess)

        theta_fit = theta_bins[:-1]
        fitted_amplitude = bisinusoidal_func(theta_fit, *popt)
        max_indices = np.argsort(fitted_amplitude)[-2:]

        direction_1 = theta_fit[max_indices[0]] * 180 / np.pi
        direction_2 = theta_fit[max_indices[1]] * 180 / np.pi
        amplitude_1 = np.max([fitted_amplitude[max_indices[0]], fitted_amplitude[max_indices[1]]])

        return amplitude_1, [direction_1, direction_2]
    else:
        return 0, []

def find_best_wavenumber_match(frames, frame_index, propagation_direction, square_size=256, overlap=0.5, strength=0.5, kl_cutoff_inf=0.01, kl_cutoff_sup=0.1, wavenumber_step=0.01, c_cutoff_inf=0, c_cutoff_sup=500, Filter=True):
    data_shape = frames.shape

    # Step 1: Pre-process the Frames
    if Filter:
        frames = filter_frames(frames, kl_cutoff_inf, kl_cutoff_sup, c_cutoff_inf, c_cutoff_sup)

    # Apply masking and get the blurred frame
    step_size = int(square_size * (1 - overlap))
    y, x = np.ogrid[:data_shape[0], :data_shape[1]]
    center_y, center_x = data_shape[0] // 2, data_shape[1] // 2
    distance_from_center = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    max_distance = strength * np.sqrt((data_shape[0] / 2)**2 + (data_shape[1] / 2)**2)
    gradient_mask = 1 - np.clip(distance_from_center / max_distance, 0, 1)
    frame = frames[:, :, frame_index] * gradient_mask

    k = fftshift(np.fft.fftfreq(data_shape[1], d=1/data_shape[1]))
    l = fftshift(np.fft.fftfreq(data_shape[0], d=1/data_shape[0]))
    k, l = np.meshgrid(k, l)
    k = k / data_shape[1]
    l = l / data_shape[0]
    radius = np.sqrt(k**2 + l**2)
    theta = np.arctan2(l, k)

    blurred_frame_fft = fftshift(fftn(frame))
    amplitude = np.abs(blurred_frame_fft)

    # Step 2: Initial Wavenumber Analysis
    wavenumber_values = np.arange(kl_cutoff_inf, kl_cutoff_sup, wavenumber_step)
    best_wavenumber = None
    best_amplitude = -np.inf
    closest_direction_diff = np.inf

    for wavenumber in wavenumber_values:
        mean_amplitude, directions = analyze_amplitude_directions(amplitude, radius, theta, wavenumber, wavenumber_step)
        for direction in directions:
            direction_diff = abs(direction - propagation_direction) % 360
            if direction_diff < closest_direction_diff:
                closest_direction_diff = direction_diff
                best_wavenumber = wavenumber
                best_amplitude = mean_amplitude

    # Step 3: Iterative Refinement
    refinement_steps = 0
    while wavenumber_step > 0.001:  # Adjust the threshold as needed
        refinement_steps += 1
        kl_cutoff_inf = max(best_wavenumber - wavenumber_step, kl_cutoff_inf)
        kl_cutoff_sup = min(best_wavenumber + wavenumber_step, kl_cutoff_sup)
        wavenumber_values = np.arange(kl_cutoff_inf, kl_cutoff_sup, wavenumber_step)
        wavenumber_step /= 2

        for wavenumber in wavenumber_values:
            mean_amplitude, directions = analyze_amplitude_directions(amplitude, radius, theta, wavenumber, wavenumber_step)
            for direction in directions:
                direction_diff = abs(direction - propagation_direction) % 360
                if direction_diff < closest_direction_diff:
                    closest_direction_diff = direction_diff
                    best_wavenumber = wavenumber
                    best_amplitude = mean_amplitude

    # Step 4: Final Wavenumber Match
    return best_wavenumber, best_amplitude, refinement_steps

def visualize_results_SFilt_Best_Wavenumber(frames, frame_index, propagation_direction, best_wavenumber, best_amplitude, avg_shift, 
                      square_number, square_size=256, overlap=0.5, strength=0.5,
                       kl_cutoff_inf=0.01, kl_cutoff_sup=0.1, c_cutoff_inf=0, 
                       c_cutoff_sup=500, Filter=True,Amp=False, log=True, Km=True):
    
    data_shape = frames.shape

    # Filter the frame first
    if Filter:
        # Perform FFT on the entire frame
        fft_data = np.fft.fftn(frames)
        fft_data = np.fft.fftshift(fft_data)

        # Calculate kl unitlessly for each pixel
        k = np.fft.fftshift(np.fft.fftfreq(data_shape[1], d=1/data_shape[1])) 
        l = np.fft.fftshift(np.fft.fftfreq(data_shape[0], d=1/data_shape[0]))
        k, l = np.meshgrid(k, l)

        k = k/(data_shape[1])
        l = l/(data_shape[0])


        kl = np.sqrt(k**2 + l**2)
        kl = np.repeat(kl[:, :,np.newaxis], data_shape[2], axis=2)


        # Frequency components as unitless (normalized index positions)
        f = np.fft.fftfreq(data_shape[2])
        f = np.repeat(f[np.newaxis, :], data_shape[1], axis=0)
        f = np.repeat(f[np.newaxis, :, :], data_shape[0], axis=0)
        f = np.fft.fftshift(f)

        # Calculate unitless phase speed
        kl[kl == 0] = np.inf  # Avoid division by zero
        c = np.abs(f) / kl  # c is unitless
        kl[kl == np.inf] = 0


        filter_mask = (c >= c_cutoff_inf) & (c <= c_cutoff_sup) & (kl >= kl_cutoff_inf) & (kl <= kl_cutoff_sup) 

        # Apply the mask to the FFT data

        fft_data[~filter_mask] = 0

        fft_data = np.fft.ifftshift(fft_data)
        fft_data = np.fft.ifftn(fft_data)
        frames = np.real(fft_data)

        if  Amp == True:
            frames = np.abs(frames)
        
    frame = frames[:,:,frame_index]
    data_shape = frame.shape

    step_size = int(square_size * (1 - overlap))
    num_squares_per_row = (data_shape[1] - square_size) // step_size + 1
    
    i = (square_number // num_squares_per_row) * step_size
    j = (square_number % num_squares_per_row) * step_size

    # Create a mask for the entire frame
    mask = np.zeros(data_shape)

    # Create a 2D radial gradient mask for the area outside the square
    y, x = np.ogrid[:data_shape[0], :data_shape[1]]
    center_y, center_x = i + square_size // 2, j + square_size // 2
    distance_from_center = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    max_distance = strength * np.sqrt((data_shape[0] / 2)**2 + (data_shape[1] / 2)**2)
    gradient_mask = 1 - np.clip(distance_from_center / max_distance, 0, 1)

    # Apply the gradient mask to the area outside the central square
    mask[i:i+square_size, j:j+square_size] = 1
    mask = np.maximum(mask, gradient_mask)

    # Apply the mask to the frame
    blurred_frame = frame * mask

    # Plot the blurred frame and the corresponding rose plot
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    ax_blur = axs[0]
    ax_polar = fig.add_subplot(122, polar=True)  # Ensure the second subplot is a polar plot
    cmap = plt.get_cmap('coolwarm')

    direction = propagation_direction, 
    amplitude = best_amplitude
    wavenumber= best_wavenumber 

    if log:
        # Apply log scale to amplitudes
        amplitude = np.log(amplitude + 1)  # Adding 1 to avoid log(0)

    # Wrap angles to the range [0, 360] degrees
    direction_wrapped = np.mod(direction, 360)

    # Plot the histogram for both directions with the same amplitude
    for d in direction_wrapped:
        if Km == True:
            ax_polar.bar(np.deg2rad(d), amplitude, width=np.deg2rad(10), color='red' , alpha=0.6, edgecolor='k', label=f'{2*1/wavenumber:.3f}' if d == direction_wrapped[0] else "")
        else:
            ax_polar.bar(np.deg2rad(d), amplitude, width=np.deg2rad(10), color='red' , alpha=0.6, edgecolor='k', label=f'{2*1/wavenumber:.3f}' if d == direction_wrapped[0] else "")

    ax_polar.set_title('Circular Spectrum of Directional Angles')
    ax_polar.set_theta_zero_location('N')
    ax_polar.set_theta_direction(-1)
    ax_polar.legend(loc='upper right', bbox_to_anchor=(1.1, 1.1))

    ax_blur.pcolormesh(blurred_frame, cmap='gray')
    ax_blur.set_title(f'Square {square_number} Blurred Frame')

    # Calculate the arrow starting point and length
    i = (square_number // num_squares_per_row) * step_size
    j = (square_number % num_squares_per_row) * step_size    
    start_x = j + square_size // 2
    start_y = i + square_size // 2
    end_x = start_x + avg_shift[1] * 0.2  # Scale factor for visualization
    end_y = start_y + avg_shift[0] * 0.2  # Scale factor for visualization

    # Plot the arrow
    ax_blur.arrow(start_x, start_y, end_x - start_x, end_y - start_y, color='red', head_width=20, head_length=30, linewidth=2)

    plt.show()