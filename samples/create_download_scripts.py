with open('healthy_beta.txt', 'r', encoding='utf-8') as f:
    healthy_beta = [line.strip() for line in f.readlines()]

with open('t2d_beta.txt', 'r', encoding='utf-8') as f:
    t2d_beta = [line.strip() for line in f.readlines()]

filepath = "shell_scripts/ena_manifest_PRJEB15401.sh"
with open(filepath, 'r', encoding='utf-8') as f:
    full_download_list = [line.strip() for line in f.readlines()]

with open("download_healthy_beta_samples.sh", "w", encoding="utf-8") as f:
    for line in full_download_list:
        _, filename = line.rsplit("/", 1)
        sample_name, _ = filename.split(".", 1)
        if sample_name in healthy_beta:
            f.write(line + "\n")

with open("download_t2d_beta_samples.sh", "w", encoding="utf-8") as f:
    for line in full_download_list:
        _, filename = line.rsplit("/", 1)
        sample_name, _ = filename.split(".", 1)
        if sample_name in t2d_beta:
            f.write(line + "\n")


print("Success")








