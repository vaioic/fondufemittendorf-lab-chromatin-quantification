from shared import tem_chromatin_analyzer as TCA

# TCA.process_image(
#     "../data/Sample #1/1.jpg",
#     "../QuPath/Sample 1/export/1-labels.png",
#     "../processed/2026-09-24 Dev",
# )

TCA.process_folder("../QuPath/Sample 1", "../processed/2026-09-25 Dev/Sample 1")
TCA.process_folder("../QuPath/Sample 2", "../processed/2026-09-25 Dev/Sample 2")
