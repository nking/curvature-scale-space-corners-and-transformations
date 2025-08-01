import torch
import torchvision
from torchvision.models.detection.faster_rcnn import FastRCNNPredictor
from torchvision.models.detection.mask_rcnn import MaskRCNNPredictor
from torchvision.transforms import v2 as T
import os
import json

def get_faster_rcnn(num_classes_custom) -> torch.nn.Module:
    # load a model pre-trained on COCO
    #Faster R-CNN predicts both bounding boxes and class scores for potential objects in the image
    model = torchvision.models.detection.fasterrcnn_resnet50_fpn_v2(weights="DEFAULT")
    
    # Freeze all the layers of the pre-trained model
    for param in model.parameters():
        param.requires_grad = False
    
    # get number of input features for the classifier
    in_features = model.roi_heads.box_predictor.cls_score.in_features
    print(f'in_features={in_features}')

    # replace the pre-trained head with a new one
    model.roi_heads.box_predictor = FastRCNNPredictor(in_features, num_classes_custom)

    return model

def get_faster_rcnn_with_masks(num_classes_custom) -> torch.nn.Module:
    # load a model pre-trained on COCO
    model = torchvision.models.detection.maskrcnn_resnet50_fpn(weights="DEFAULT")
    
    # get number of input features for the classifier
    in_features = model.roi_heads.box_predictor.cls_score.in_features
    print(f'in_features={in_features}')

    # replace the pre-trained head with a new one
    model.roi_heads.box_predictor = FastRCNNPredictor(in_features, num_classes_custom)

    in_features_mask = model.roi_heads.mask_predictor.conv5_mask.in_channels
    hidden_layer = 256
    print(f'hidden_layer={hidden_layer}')
    # and replace the mask predictor with a new one
    model.roi_heads.mask_predictor = MaskRCNNPredictor(
        in_features_mask,
        hidden_layer,
        num_classes_custom
    )

    return model

def get_transform(train:bool=True):
    transforms = []
    if train:
        transforms.append(T.RandomHorizontalFlip(0.5))
    transforms.append(T.ToDtype(torch.float, scale=True))
    transforms.append(T.ToPureTensor())
    return T.Compose(transforms)

def quantize_model(model):
    import torch.quantization
    # Apply dynamic quantization to reduce model weights to lower precision and activations are converted to
    # lower precision on the fly.  this is useful when need to reduce memory access times, such as in LSTMs
    # and Transformers with small batch sizes.
    #dynamic quantization (weights quantized with activations read/stored in floating point and quantized 
    # for compute)
    #https://docs.pytorch.org/tutorials/recipes/recipes/dynamic_quantization.html
    # to find which engines are supported: torch.backends.quantized.supported_engines
    # first try, runtime error due to quantized::linear_prepack operation missing or misconfigured
    # then saw issue and gemin advice to set this to the listed engine qnnpack
    torch.backends.quantized.engine = 'qnnpack'
    quantized_model = torch.quantization.quantize_dynamic(
        model, {torch.nn.Linear}, dtype=torch.qint8
    )
    return quantized_model

def get_json_param_filepath(model_dir, model_name):
    return os.path.join(model_dir, f"{model_name}_params.json")

def get_quantized_json_param_filepath(model_dir, model_name):
    return os.path.join(model_dir, f"{model_name}_quan_params.json")

def get_saved_model_filepath(model_dir, model_name):
    return os.path.join(model_dir, f"{model_name}.pth")

def get_quantized_saved_model_filepath(model_dir, model_name):
    return os.path.join(model_dir, f"{model_name}_quan.pth")
    
def save(model, params, model_dir, model_name):
    '''
    Args:
        model:

        params: containing parameters like num_classes_custom
        
        model_dir: e.g. saved_models

        model_name: e.g. churn_model_v20250801_01
    '''
    #params = convert_float32_to_float64(params)
    with open(get_json_param_filepath(model_dir, model_name), 'w') as f:
        json.dump(params, f, indent=4)
    torch.save(model.state_dict(), get_saved_model_filepath(model_dir, model_name))

def _save_quantized(model, params, model_dir, model_name):
    #params = convert_float32_to_float64(params)
    params = params.copy()
    if 'optimizer' in params:
        params.pop('optimizer')
    if 'criterion' in params:
        params.pop('criterion')
    with open(get_quantized_json_param_filepath(model_dir, model_name), 'w') as f:
        json.dump(params, f, indent=4)
    q_model = quantize_model(model)
    torch.save(q_model.state_dict(), get_quantized_saved_model_filepath(model_dir, model_name))

def load_model(model_dir, model_name, device):
    '''
    Args:        
        model_dir: e.g. saved_models

        model_name: e.g. churn_model_v20250801_01

        device
    '''
    with open(get_json_param_filepath(model_dir, model_name), 'r') as f:
        params = json.load(f)

    num_classes_custom = params['num_classes_custom']
    
    model = get_faster_rcnn(num_classes_custom)
    
    checkpoint = torch.load(get_saved_model_filepath(model_dir, model_name), \
        map_location=device) # Use 'cuda' if loading to GPU

    model.load_state_dict(checkpoint)

    model.eval()
    
    return model, params
        
def _load_quantized_model(model_dir, model_name):
    with open(get_quantized_json_param_filepath(model_dir, model_name), 'r') as f:
        params = json.load(f)
    num_classes_custom = params['num_classes_custom']

    model = get_faster_rcnn(num_classes_custom)
    model.qconfig = torch.quantization.get_default_qconfig('qnnpack')
    #model_prepared = torch.quantization.prepare(model)
    model_prepared = torch.quantization.prepare_qat(model)
    model_quantized = torch.quantization.convert(model_prepared)
    model_quantized.load_state_dict(torch.load(get_quantized_saved_model_filepath(model_dir, model_name)))
    return model, params

class EarlyStopping:
    '''
    from https://medium.com/biased-algorithms/a-practical-guide-to-implementing-early-stopping-in-pytorch-for-model-training-99a7cbd46e9d
    '''
    def __init__(self, patience=5, delta=0, verbose=False):
        self.patience = patience
        self.delta = delta
        self.verbose = verbose
        self.best_loss = None
        self.no_improvement_count = 0
        self.stop_training = False
    
    def check_early_stop(self, val_loss):
        if self.best_loss is None or val_loss < self.best_loss - self.delta:
            self.best_loss = val_loss
            self.no_improvement_count = 0
        else:
            self.no_improvement_count += 1
            if self.no_improvement_count >= self.patience:
                self.stop_training = True
                if self.verbose:
                    print("Stopping early as no improvement has been observed.")
                    
