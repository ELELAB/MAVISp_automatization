import argparse
import requests
import csv

def fetch_ted_domains(uniprot_id):
    url = f"https://ted.cathdb.info/api/v1/uniprot/summary/{uniprot_id}?skip=0&limit=100"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        return response.json().get("data", [])
    except requests.RequestException as e:
        print(f"Error fetching data for UniProt ID {uniprot_id}: {e}")
        return []

def find_matching_domains(uniprot_id):
    domains = fetch_ted_domains(uniprot_id)
    matched = []

    for domain in domains:
        chop_str = domain.get("chopping", "")
        if not chop_str:
            continue

        cath_label_raw = domain.get("cath_label", "").strip()
        cath_label = "" if cath_label_raw == "-" or not cath_label_raw else cath_label_raw

        matched.append({
            "TED_id": domain.get("ted_id"),
            "TED_boundaries": chop_str,
            "CATH_label": cath_label
        })
    return matched

def main():
    parser = argparse.ArgumentParser(description="Query TED domain annotations for a UniProt ID.")
    parser.add_argument("-id", "--uniprot_id", help="UniProt ID (e.g., P36896)", required=True)
    parser.add_argument("-o", "--output", help="Output CSV file path", default="results.csv")
    args = parser.parse_args()

    matching_domains = find_matching_domains(args.uniprot_id)

    with open(args.output, mode="w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["TED_id", "TED_boundaries", "CATH_label"])
        writer.writeheader()
        writer.writerows(matching_domains)

if __name__ == "__main__":
    main()
