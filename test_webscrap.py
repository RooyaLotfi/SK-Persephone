import os
from Utils import pull_html_data


def main():
    directory = os.path.join('data', 'HLS')
    archive_url = "https://192.168.1.12/data-04/HLS/"

    pull_html_data.SATELLITE = {'L30', 'S30'}
    # First three characters of the tiles
    pull_html_data.TILE = {'14ULA'}
    # Extension of files to be downloaded
    pull_html_data.EXTENSION = {'.tif'}
    pull_html_data.YEAR = {'2019'}
    # The days should fall between day 150 and 180 of the year
    pull_html_data.DAY_INTERVAL = [150, 180]
    pull_html_data.USER_NAME = 'airm-user'
    pull_html_data.PASSWORD = 'bh7S-p0Ml-Wsrt-Ler3'

    pull_html_data.crawl_web(archive_url, directory)


if __name__ == "__main__":
    main()
