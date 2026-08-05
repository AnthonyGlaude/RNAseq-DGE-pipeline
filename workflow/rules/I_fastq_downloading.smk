import pandas as pd
from pathlib import Path


rule fastq_download:
    input:
        urls=config["fastq"]["url_list"]
    output:
        done=touch(config["fastq"]["download_done"])
    params:
        fastq_dir=config["path"]["fastq"]
    log: 
        "logs/fastq_download.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.fastq_dir}
        mkdir -p "$(dirname {log})"

        wget \
            --continue \
            --input-file {input.urls} \
            --directory-prefix {params.fastq_dir} \
            > {log} 2>&1

        if ! find {params.fastq_dir} \
            -maxdepth 1 \
            -name "*.fastq.gz" \
            -print \
            -quit \
            | grep -q .; then

            echo "Aucun FASTQ trouvé dans {params.fastq_dir}" >> {log}
            exit 1
        fi
        """


rule rename_fastq:
    input:
        rename_table=config["fastq"]["rename_tsv"],
        download_done=config["fastq"]["download_done"]
    output:
        done=touch(config["fastq"]["rename_done"])
    params:
        fastq_dir=config["path"]["fastq"]
    log: 
        "logs/rename_fasta.log"
    run:
        fastq_dir = Path(params.fastq_dir)
        Path(log[0]).parent.mkdir(parents=True, exist_ok=True)

        df = pd.read_csv(
            input.rename_table,
            sep=";",
            dtype=str
        )

        required_columns = {
            "sample",
            "condition",
            "genotype"
        }

        missing_columns = required_columns - set(df.columns)

        if missing_columns:
            raise ValueError(
                "Colonnes manquantes : "
                + ", ".join(sorted(missing_columns))
            )

        for column in required_columns:
            df[column] = df[column].str.strip()

        errors = []

        with open(log[0], "w") as logfile:
            for row in df.itertuples(index=False):
                for read in ("R1", "R2"):
                    src = fastq_dir / (
                        f"{row.sample}_{read}.fastq.gz"
                    )

                    dst = fastq_dir / (
                        f"{row.condition}_{row.genotype}_"
                        f"{row.sample}_{read}.fastq.gz"
                    )

                    if src.exists():
                        if dst.exists():
                            errors.append(
                                f"Destination déjà existante : {dst}"
                            )
                            continue

                        src.rename(dst)
                        logfile.write(
                            f"{src.name} -> {dst.name}\n"
                        )

                    elif dst.exists():
                        logfile.write(
                            f"Déjà renommé : {dst.name}\n"
                        )

                    else:
                        errors.append(
                            f"Fichier absent : {src}"
                        )

        if errors:
            raise RuntimeError(
                "Erreur pendant le renommage :\n"
                + "\n".join(errors)
            )
