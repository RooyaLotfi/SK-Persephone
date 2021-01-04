import os

import requests
# turn off warnings for the unverified SSL connection
import urllib3
from bs4 import BeautifulSoup


def get_links(archive_url):
    # create response object
    r = requests.get(archive_url, auth=('airm-user', 'bh7S-p0Ml-Wsrt-Ler3'), verify=False)

    # create beautiful-soup object
    soup = BeautifulSoup(r.content, 'html.parser')

    # find all links on web-page
    links = soup.findAll('a')

    # filter the link sending with .mp4
    links = [archive_url + link['href'] for link in links if (link['href'].endswith('.hdf') or link['href'].endswith('.hdr'))]

    return links


def download_series(links):
    for link in links:

        '''iterate through all links in links  
        and download them one by one'''

        # obtain filename by splitting url and getting
        # last string
        file_name = link.split('/')[-1]

        print("Downloading file:%s" % file_name)

        # create response object
        r = requests.get(link, auth=('airm-user', 'bh7S-p0Ml-Wsrt-Ler3'), verify=False)

        # download started
        with open(r'/Users/roya.lotfi/AIRM_Project/SK-Persephone/data/hls_tmp/'+file_name, 'wb') as f:
            f.write(r.content)

        print("%s downloaded!\n" % file_name)

    print("All tif files downloaded!")


def download_all(archive_url, directory):
    urllib3.disable_warnings()
    if not os.path.exists(directory):
        os.makedirs(directory)
    links = get_links(archive_url)
    download_series(links)


def main():
    directory = 'hls_tmp/'
    archive_url = "https://192.168.1.12/data-04/hls_tmp/12UWA/2019/L30/"
    download_all(archive_url, directory)


if __name__ == "__main__":
    main()
