import os
import pathlib
import matplotlib
import matplotlib.pyplot as plt
import io
import imageio
import scipy.misc
import numpy as np
from six import BytesIO
import PIL
from PIL import Image, ImageDraw, ImageFont
from IPython.display import display, Javascript
from IPython.display import Image as IPyImage
import torch
from torch.utils.data import Dataset
from torchvision import tv_tensors
from torchvision.tv_tensors import BoundingBoxFormat
from torchvision.utils import draw_bounding_boxes, draw_segmentation_masks
from torchvision.io import read_image
from torchvision.ops.boxes import masks_to_boxes
from torchvision.transforms.v2 import functional as F


#import tensorflow as tf

FIGSIZE = (8, 6)
THRESH = 0.5

#NOTE:  PyTorch labels use class 0 as background. If your dataset does not contain the 
# background class, you should not have 0 in your labels
# ==> everything below here gets +1 added to its labels

class_mapping = {0:'cupcake', 1:'euclair', 2:'icecream', 3:'gingerbread_man',
  4:'icecream_sandwich', 5:'honey_comb', 6:'kitkat',
  7:'jellybean', 8:'donut'}

class_mapping_1 = {1:'cupcake', 2:'euclair', 3:'icecream', 4:'gingerbread_man',
  5:'icecream_sandwich', 6:'honey_comb', 7:'kitkat',
  8:'jellybean', 9:'donut'}

#"xyxy"	Boxes defined by (x_min, y_min, x_max, y_max) — corners of the box
#"yxyx"	Boxes defined by (y_min, x_min, y_max, x_max) — same as above, but reordered
#"xywh"	Boxes defined by (x_center, y_center, width, height) — center + dimensions
#"center_xywh"	Used internally for encoding offsets during training

def make_bb(y_min, x_min, y_max, x_max, w, h) :
  '''
  for pytorch
  "xyxy"	Boxes defined by (x_min, y_min, x_max, y_max)
  '''
  return [x_min, y_min, x_max, y_max]

def load_image_into_numpy_array(path : str):
  """Load an image from path into a numpy array.

  Args:
    path: a file path.

  Returns:
    uint8 numpy array with shape (img_height, img_width, 3)
  """
  img_data = tf.io.gfile.GFile(path, 'rb').read()
  image = Image.open(BytesIO(img_data))
  (im_width, im_height) = image.size
  return np.array(image.getdata()).reshape((im_height, im_width, 3)).astype(np.uint8)

class PyTorchDataset(torch.utils.data.Dataset):
    def __init__(self, image_filepaths, boxes, class_labels : list, image_ids, transforms):
        if len(image_filepaths) != len(boxes) or len(image_filepaths) != len(class_labels):
            raise ValueError("image_filepaths and labels and class_labels have the same length.")
        self.image_filepaths = image_filepaths
        self.boxes = boxes
        self.class_labels = class_labels
        self.image_ids = image_ids
        self.transforms = transforms
        
    def __len__(self):
        return len(self.image_filepaths)
        
    def __getitem__(self, idx):
        ''' return tuple of
        image: 
            torchvision.tv_tensors.Image of shape [3, H, W], a pure tensor, 
            or a PIL Image of size (H, W)
        target: 
            a dict containing the following fields
              boxes:
                torchvision.tv_tensors.BoundingBoxes of shape [N, 4]
                the coordinates of the N bounding boxes in [x0, y0, x1, y1] format, 
                ranging from 0 to W and 0 to H
              labels:
                integer torch.Tensor of shape [N]: the label for each bounding box. 
                0 represents always the background class.
              image_id, 
                int: an image identifier. It should be unique between all the images in the dataset, 
                and is used during evaluation
              area:
                float torch.Tensor of shape [N]: the area of the bounding box. 
                This is used during evaluation with the COCO metric, to separate the metric scores 
                between small, medium and large boxes.
              iscrowd:
                uint8 torch.Tensor of shape [N]: instances with iscrowd=True will be ignored 
                during evaluation.
              masks:
                  (optional)
                  torchvision.tv_tensors.Mask of shape [N, H, W]: the segmentation masks for each one 
                  of the objects
        '''
        image = tv_tensors.Image(read_image(self.image_filepaths[idx]))
        bounds = self.boxes[idx]
        labels = torch.tensor(self.class_labels[idx])
        image_id = torch.tensor(self.image_ids[idx])
        area = (bounds[:, 3] - bounds[:, 1])*(bounds[:, 2] - bounds[:, 0])
        iscrowd = torch.tensor([0 for i in range(len(bounds))])
        target = {'boxes':bounds, 'labels':labels, 'image_id' : image_id, 'area':area, \
            'iscrowd':iscrowd}
        if self.transforms is not None:
            image, target = self.transforms(image, target)
        return image, target
        
class_ids_01 = [0, 1, 2, 3, 4, 5]
bb_01 = [make_bb(134,     18,      191,      70, 640, 279),
    make_bb(176,    101,      244,     229, 640, 279),
    make_bb(129,    125,      188,     175, 640, 279),
    make_bb(122,    212,      193,     262, 640, 279),
    make_bb(101,   280,      213,     365, 640, 279),
    make_bb(128,   391,      190,     455, 640, 279),
]

class_ids_02 = [0, 1, 2, 3]
bb_02 = [make_bb(125,   128,        213,    212, 640, 427),
    make_bb(186,   61,         358,    207, 640, 427),
    make_bb(105,  274,         233,    376, 640, 427),
    make_bb(59,   493,         272,    594, 640, 427),
]

class_ids_03 = [1, 3, 4, 5, 6, 7]
bb_03 = [make_bb(396,   278,       467,     497, 1280, 960),
    make_bb(234,    57,       525,     303, 1280, 960),
    make_bb(227,  501,       494,     740, 1280, 960),
    make_bb(279,   702,       447,     788, 1280, 960),
    make_bb(228,   494,       590,     714, 1280, 960),
    make_bb(255,   932,       477,    1046, 1280, 960),
]

class_ids_04 = [0,1,2,3]
bb_04 = [make_bb(309,   311,       615,    605, 1280, 960),
    make_bb(456,   998,       507,   1131, 1280, 960),
    make_bb(329,   671,       575,    848, 1280, 960),
    make_bb(318,   875,       552,   1041, 1280, 960),
]

class_ids_05 = [3, 4]
bb_05 = [make_bb(210,   107,       413,   177, 450, 600),
    make_bb(161,   152,      507,   480, 450, 600),
]

class_ids_06 = [0, 6, 8]
bb_06 = [make_bb(33,     0,       121,    41, 320, 181),
    make_bb(44,   239,        97,   261, 320, 181),
    make_bb(0,      58,      180,   244, 320, 181)
]

class_ids_07 = [7, 8]
bb_07 = [make_bb(15,    151,       197,    274, 320, 213),
    make_bb(64,     80,       110,    112, 320, 213)
]

def read_android_image(idx):
    img_path = os.path.join("../testresources", 'android_statues_0' + str(idx) + '.jpg')
    img = read_image(img_path)
    return img

def load_statues(transforms, add_one_to_classes:bool=True):
    train_files = [os.path.join("../testresources", 'android_statues_0' + str(i) + '.jpg') \
        for i in range(1, 8)]
    #Expected target boxes to be a tensor of shape [N, 4], got torch.Size([1, 6, 4]).
    bboxes = []
    bboxes.append(tv_tensors.BoundingBoxes(torch.tensor(bb_01, dtype=torch.float32),\
        format=BoundingBoxFormat.XYXY, canvas_size=[279, 640]))
    bboxes.append(tv_tensors.BoundingBoxes(torch.tensor(bb_02, dtype=torch.float32),\
        format=BoundingBoxFormat.XYXY, canvas_size=[427, 640]))
    bboxes.append(tv_tensors.BoundingBoxes(torch.tensor(bb_03, dtype=torch.float32),\
        format=BoundingBoxFormat.XYXY, canvas_size=[960, 1280]))
    bboxes.append(tv_tensors.BoundingBoxes(torch.tensor(bb_04, dtype=torch.float32),\
        format=BoundingBoxFormat.XYXY, canvas_size=[960, 1280]))
    bboxes.append(tv_tensors.BoundingBoxes(torch.tensor(bb_05, dtype=torch.float32),\
        format=BoundingBoxFormat.XYXY, canvas_size=[600, 450]))
    bboxes.append(tv_tensors.BoundingBoxes(torch.tensor(bb_06, dtype=torch.float32),\
        format=BoundingBoxFormat.XYXY, canvas_size=[181, 320]))
    bboxes.append(tv_tensors.BoundingBoxes(torch.tensor(bb_07, dtype=torch.float32),\
        format=BoundingBoxFormat.XYXY, canvas_size=[213, 320]))

    if add_one_to_classes:
        class_ids = [(np.array(class_ids_01)+1).tolist(), (np.array(class_ids_02)+1).tolist(),\
            (np.array(class_ids_03)+1).tolist(), (np.array(class_ids_04)+1).tolist(),\
            (np.array(class_ids_05)+1).tolist(), (np.array(class_ids_06)+1).tolist(),\
            (np.array(class_ids_07)+1).tolist()]
        category_index = {i: {'id':i, 'name':cls} for i, cls in class_mapping_1.items()}
        num_classes = len(category_index) + 1
    else:
        class_ids = [class_ids_01, class_ids_02, class_ids_03, class_ids_04,
            class_ids_05, class_ids_06, class_ids_07]
        category_index = {i: {'id':i, 'name':cls} for i, cls in class_mapping.items()}
        num_classes = len(category_index)
    image_ids = [1,2,3,4,5,6,7] 
    ds = PyTorchDataset(train_files, bboxes, class_ids, image_ids, transforms)
    return category_index, ds, None

def load_gbman_3(transforms, add_one_to_classes:bool=True):
    train_files = [os.path.join("../testresources", 'android_statues_0' + str(i) + '.jpg') \
        for i in [1,2,4]]
    test_files = [os.path.join("../testresources", 'android_statues_0' + str(i) + '.jpg') \
        for i in [3, 5]]
    #error here Tensor[N, 4]
    #Expected target boxes to be a tensor of shape [N, 4], got torch.Size([1, 1, 4]
    bboxes = []
    bboxes.append(tv_tensors.BoundingBoxes(torch.tensor(bb_01[3], dtype=torch.float32),\
        format=BoundingBoxFormat.XYXY, canvas_size=[279, 640]))
    bboxes.append(tv_tensors.BoundingBoxes(torch.tensor(bb_02[3], dtype=torch.float32),\
        format=BoundingBoxFormat.XYXY, canvas_size=[427, 640]))
    bboxes.append(tv_tensors.BoundingBoxes(torch.tensor(bb_04[3], dtype=torch.float32),\
        format=BoundingBoxFormat.XYXY, canvas_size=[960, 1280]))
    class_ids = [[0], [0], [0]]
    image_ids = [1,2,4]
              
    bboxes_test=[]
    bboxes_test.append(tv_tensors.BoundingBoxes(torch.tensor(bb_03[1], dtype=torch.float32),\
        format=BoundingBoxFormat.XYXY, canvas_size=[960, 1280]))
    bboxes_test.append(tv_tensors.BoundingBoxes(torch.tensor(bb_05[0], dtype=torch.float32),\
        format=BoundingBoxFormat.XYXY, canvas_size=[600, 450]))
    class_ids_test = [[0], [0]]
    image_ids_test = [3,5]

    if add_one_to_classes:
        class_ids = [[1], [1], [1]]
        class_ids_test = [[1], [1]]
        category_index = {1: {'id': 1, 'name': 'gingerbread_man'}}
        num_classes = len(category_index) + 1
    else:
        category_index = {0: {'id': 0, 'name': 'gingerbread_man'}}
        num_classes = len(category_index)
    
    ds = PyTorchDataset(train_files, bboxes, class_ids, image_ids, transforms)
    ds_test = PyTorchDataset(test_files, bboxes_test, class_ids_test, image_ids_test, transforms)
    return category_index, ds, ds_test
       
def get_video_snapshot_filepaths() -> []:
    #!pip install -U --pre "yt-dlp[default]"
    #!pip install opencv-python
    #!pip install --upgrade yt-dlp
    #!apt -y install ffmpeg lame  #or install w/ brew or other package manager
    #!pip install imutils
    import imutils
    import cv2
    
    video_url = "https://www.youtube.com/watch?v=BRKLw_16Lac"
    
    data_dir = os.path.join(os.getcwd(), '..', 'bin')
    
    def get_frame(time, frame_count, filepath) :
        vid_cap.set(cv2.CAP_PROP_POS_MSEC, time)
        # youtube frame rate options: 24 to 60 frames/sec
        frame_det, frame = vid_cap.read()
        print(f'read correctly={frame_det}')
        if frame_det:
            frame = imutils.resize(frame, width=640) #smaller width of 128?
            print(f'writing to {filepath}')
            return cv2.imwrite(filepath, frame) #return True if written
        else:
            return None
    
    stream_uri = os.path.join(data_dir, "android_statues.mp4")
    
    #https://pypi.org/project/yt-dlp/
    #Video Format Options
    # --check-formats
    
    #!yt-dlp $video_url --list-formats
    
    if not os.path.exists(stream_uri):
        import subprocess
        print(f'downloading youtube file video_url={video_url}\nstream_uri={stream_uri}')
        # default is ffmpeg,  can choose mp4 with -f
        #!yt-dlp $video_url -vU -f mp4 -o $stream_uri
        os.environ["video_url"] = video_url
        os.environ["stream_uri"] = stream_uri
        shell_path = os.getenv('SHELL', '/bin/bash')
        process = subprocess.Popen([shell_path, '-i', '-l'], \
            stdin=subprocess.PIPE, \
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate(input='yt-dlp $video_url -vU -f mp4 -o $stream_uri') 
        print(":", stdout.strip())
        #subprocess.run(['/bin/bash', 'yt-dlp', video_url, '-vU', '-f', 'mp4', '-o', stream_uri])
    else:
        print(f'have already downloaded YouTube mp4 file')

    filepaths = []
    
    data_dir_snaps = os.path.join(data_dir, "snapshots")
    if os.path.exists(data_dir_snaps) and os.listdir(data_dir_snaps):
        for f in os.listdir(data_dir_snaps):
            chk = os.path.join(data_dir_snaps, f)
            if os.path.isfile(chk):
                filepaths.append(chk)
    else:
        print(f'extract video frames')
        if not os.path.exists(data_dir_snaps):
            os.mkdir(data_dir_snaps)
        # the video is 35 sec long
        vid_cap = cv2.VideoCapture(stream_uri)
        if not vid_cap.isOpened():
            print(f'video capture is not opened yet.  try now:')
            vid_cap.open(video_url)
        print(f'vid_cap.isOpened()={vid_cap.isOpened()}')
    
        #vid_cap.set(cv2.CV_CAP_PROP_FRAME_WIDTH, 640)
    
        fps = vid_cap.get(cv2.CAP_PROP_FPS)
        print(f'fps={fps}')
        print(f'CAP_PROP_FRAME_COUNT={vid_cap.get(cv2.CAP_PROP_FRAME_COUNT)}')
        num_frames = int(vid_cap.get(cv2.CAP_PROP_FRAME_COUNT))
        if (fps > 0):
            end_time = int(num_frames/fps)*1000
        else:
            end_time = 35*1000
        print(f'end_time={end_time}')

        start_time = 0
        frame_rate = 2000 #in millisec = 2 sec
        frame_count = 1
        for time in range(start_time, end_time, frame_rate):
            time = round(time, 2)
            filepath = os.path.join(data_dir_snaps, str(frame_count) + '.jpg')
            result = get_frame(time, frame_count, filepath)
            frame_count += 1
            if result is not None:
                filepaths.append(filepath)
  
    return filepaths
    
    '''
    android statues 01 is from:
    https://www.flickr.com/photos/67287915@N00/8570385915
    android statues 02 is from:
    https://www.flickr.com/photos/quinnanya/5847206255
    android statues 03 and 04 are from:
    https://github.com/nking/curvature-scale-space-corners-and-transformations.git
    android statues 05 is from:
    https://commons.wikimedia.org/wiki/File:IceCream_Sandwich_%287791561448%29.jpg
    android statues 06 is from:
    https://upload.wikimedia.org/wikipedia/commons/thumb/3/3c/Sculpture_for_Android_Donut_at_Google_Mountain_View.jpg/320px-Sculpture_for_Android_Donut_at_Google_Mountain_View.jpg
    android statues 07 is from:
    https://upload.wikimedia.org/wikipedia/commons/thumb/e/ef/Android_Jelly_Bean_Lawn_Statue_%2812757851595%29.jpg/320px-Android_Jelly_Bean_Lawn_Statue_%2812757851595%29.jpg

                      ytop    xleft   ybottom  xright    w    h
      android 1 (tr)   w=640  h=279
    0 cupcake          134,     18,      191,      70,
    1 euclair          176,    101,      244,     229,
    2 icecream         129,    125,      188,     175,
    3 gingerbread_man  122,    212,      193,     262,
    4 icecream_sandwich 101,   280,      213,     365,
    5 honey_comb        128,   391,      190,     455,

      android 2 (tr)   w=640, h=427
    0 cupcake          125,   128,        213,    212
    1 euclair          186,   61,         358,    207,
    2 icecream         105,  274,         233,    376
    3 gingerbread_man  59,   493,         272,    594

      android 3 (te)  w=1280, h=960
    1 euclair          396,   278,       467,     497,
    3 gingerbread_man  234,    57,       525,     303,
    4 icecream_sandwich 227,  501,       494,     740,
            (obscured)
    5 honeycomb        279,   702,       447,     788,
    6 kitkat           228,   494,       590,     714,
    7 jellybean        255,   932,       477,    1046,

      android 4 (tr) w=1280, h=960
    0 cupcake          309,   311,       615,    605,
    1 euclair          456,   998,       507,   1131,
    2 icecream         329,   671,       575,    848,
    3 gingerbread_man  318,   875,       552,   1041,

      android 5 (te)  w=450, h=600
    3 gingerbread_man  210,   107,       413,   177,
    4 icecream_sandwich 161,   152,      507,   480,

      android 6 (tr)  w=320,  h=181
    0 cupcake           33,     0,       121,    41,
    6 kitkat            44,   239,        97,   261,
    8 donut             0,      58,      180,   244,

      android 7 (tr)  w=320,  h=213
    7 jellybean        15,    151,       197,    274,
    8 donut            64,     80,       110,    112,
                      ytop    xleft   ybottom  xright
    '''
