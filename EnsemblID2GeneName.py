import requests


def get_gene_name(transcript_id):
    ensembl_api_url = f"https://rest.ensembl.org/lookup/id/{transcript_id}?expand=1"

    response = requests.get(ensembl_api_url, headers={"Content-Type": "application/json"})

    if response.status_code == 200:
        data = response.json()
        if 'display_name' in data:
            return data['display_name']

    return None


transcript_id = "ENSG00000210195"
gene_name = get_gene_name(transcript_id)

if gene_name:
    print(f"Transcript ID: {transcript_id}")
    print(f"Gene Name: {gene_name}")
else:
    print(f"Gene name not found for Transcript ID: {transcript_id}")
