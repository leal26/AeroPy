import numpy as np
from paraview.simple import STLReader, Transform, Slice, SaveData


def cosine_spacing(start, stop, num=50, offset=0):
    # calculates the cosine spacing
    index = np.linspace(0., 1., num)
    spacing = .5*(1.-np.cos(np.pi*(index-offset)))

    points = start+spacing*(stop-start)

    return points


def get_slices(filename, directory='.\\', rotation=[0.0, -2.3067, 0.0],
               slice_coordinates=cosine_spacing(0.0001, 38.6736, 5),
               slice_direction=[1, 0, 0]):
    # create a new 'STL Reader'

    input = directory + '\\raw\\' + filename + '.stl'
    output = directory + '\\processed\\' + filename + '_'
    fuselagestl = STLReader(FileNames=[input])

    # create a new 'Transform'
    transform1 = Transform(Input=fuselagestl)

    # Properties modified on transform1.Transform
    transform1.Transform.Rotate = rotation

    for i, x in enumerate(slice_coordinates):
        # create a new 'Slice'
        slice1 = Slice(Input=transform1)
        print(vars(slice1))
        # Properties modified on slice1.SliceType
        slice1.SliceType.Origin = [x*y for y in slice_direction]

        # save data
        SaveData(output+str(i)+'.csv', proxy=slice1, Precision=12)


directory = 'D:\\GitHub\\AeroPy\\examples\\JAXA_files'
filename = 'fuselage'
get_slices(filename, directory,
           slice_coordinates=cosine_spacing(0.0001, 38.6736, 100))
