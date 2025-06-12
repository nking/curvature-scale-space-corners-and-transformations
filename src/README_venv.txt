to create a venv:
   python3 -m venv venv
   (the later can be named anything, but venv is standard)

to activate venv:
   source venv/bin/activate
   (or equiv for other than unix)

to confirm that venv is activated:
   which python

to deactivate:
   deactivate



for installing jax and jaxlib:

https://developer.apple.com/metal/jax/

python -m pip install -U pip
python -m pip install numpy wheel
python -m pip install jax-metal
python -c 'import jax; print(jax.numpy.arange(10))'\n
pip install -U jaxlib jax\n
ENABLE_PJRT_COMPATIBILITY=1 python -c 'import jax; print(jax.numpy.arange(10))'\n

