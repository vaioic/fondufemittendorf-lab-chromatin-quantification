from shared import tem_chromatin_analyzer as TCA
from shared.tem_chromatin_analyzer import combine_csv_files

TCA.process_folder("../QuPath/Sample 1", "../processed/2026-10-19/Sample 1")
TCA.process_folder("../QuPath/Sample 2", "../processed/2026-10-19/Sample 2")

TCA.process_folder("../QuPath/Sample 3", "../processed/2026-10-19/Sample 3")
TCA.process_folder("../QuPath/Sample 4", "../processed/2026-10-19/Sample 4")

TCA.process_folder("../QuPath/Sample 5", "../processed/2026-10-19/Sample 5")
TCA.process_folder("../QuPath/Sample 6", "../processed/2026-10-19/Sample 6")

combine_csv_files("../processed/2026-10-19")
