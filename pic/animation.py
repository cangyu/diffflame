import imageio

begin_idx = input('Starting index: ')
end_idx = input('Ending index: ')

image_list = ['Iteration_{}.png'.format(i) for i in range(int(begin_idx), int(end_idx)+1)]

frames = [imageio.imread(f) for f in image_list]

imageio.mimsave('diffflame.gif', frames, 'GIF', duration = 0.1)
