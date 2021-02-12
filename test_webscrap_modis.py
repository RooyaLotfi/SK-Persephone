import os
from Utils import pull_html_data


def main():
    directory = os.path.join("/Users/roya.lotfi/AIRM_Project/SK-Persephone/data/Modis_MCD43A4/2019/")
    archive_url = "https://192.168.1.12/data-04/MCD43A4/2019/"

    pull_html_data.YEAR = {'2019'}
    # The days should fall between day 150 and 180 of the year
    pull_html_data.DAY_INTERVAL = [151, 190]
    pull_html_data.TILE = {'h10v03', 'h10v04', 'h11v03', 'h11v04', 'h12v03'}
    pull_html_data.EXTENSION = {'.hdr', '.hdf', '.tif', '.xml'}
    pull_html_data.USER_NAME = 'airm-user'
    pull_html_data.PASSWORD = 'bh7S-p0Ml-Wsrt-Ler3'

    pull_html_data.crawl_web(archive_url, directory)


if __name__ == "__main__":
    main()
