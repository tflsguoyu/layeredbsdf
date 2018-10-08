import os

argv = '..\\..\\..\\dist\\mtsutil '
argv0 = 'traceray_height '
# argv1 = 'Patch_height_513.exr '
argv1 = 'scene.xml '
argv2 = '1000000 '
argv3 = 'traceray_sample.exr '
argv4 = 'traceray_eval.exr'
os.system(argv + argv0 + argv1 + argv2 + argv3 + argv4)
