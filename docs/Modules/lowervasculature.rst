================================
Lower_generation_segmentation.py
================================

Basically, this uses the same method as lung segmentation (watershed). With watershed algorithm the vessel gaps inside the lung masks are used. One step of binary opening is done on the masks to remove really small vessel gaps inside the lung masks since for our purpose we are interested in bigger and more visible vessels.
