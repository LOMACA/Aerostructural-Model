# URL Correction Function
def fix_url(relative_url):
    # Ensure the URL starts with a correct base (http or https)
    if not relative_url.startswith("http"):
        return f"https://m-selig.ae.illinois.edu/ads/{relative_url.lstrip('/')}"
    return relative_url

def download_all_airfoil_dat_files(base_url="https://m-selig.ae.illinois.edu/ads/coord_database.html",
                                   output_dir="airfoil_dat_files"):
    import os
    import time
    import requests
    from bs4 import BeautifulSoup

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    response = requests.get(base_url)
    if response.status_code != 200:
        raise Exception("Failed to retrieve airfoil page")

    soup = BeautifulSoup(response.text, 'html.parser')

    dat_links = soup.select("a[href$='.dat']")

    for dat_link in dat_links:
        airfoil_name = dat_link.text.strip()
        dat_file_url = fix_url(dat_link['href'])

        print(f"Downloading DAT file from URL: {dat_file_url}")

        response = requests.get(dat_file_url)
        if response.status_code != 200:
            print(f"Failed to download DAT file for {airfoil_name}")
            continue

        dat_file_path = os.path.join(output_dir, f"{airfoil_name}.dat")
        with open(dat_file_path, 'wb') as f:
            f.write(response.content)

        print(f"Downloaded DAT file for {airfoil_name} to {dat_file_path}")

        time.sleep(1)

download_all_airfoil_dat_files()