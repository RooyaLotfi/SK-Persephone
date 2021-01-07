"""
Designed to get the HLS data from server: https://192.168.1.12/data-04/hls_tmp/
"""

import os
import requests
# turn off warnings for the unverified SSL connection
import urllib3
from bs4 import BeautifulSoup

# Default values for satellite, tile, extension of files to download, and year
SATELLITE = {'L30', 'S30'}
# First three characters of the tiles
TILE = {'13UBS', '13UBR'}
# Extension of files to be downloaded
EXTENSION = {'.hdr', '.hdf'}
YEAR = {'2019'}
USER_NAME = 'airm-user'
PASSWORD = 'bh7S-p0Ml-Wsrt-Ler3'


def get_links(archive_url):
    """

    :param archive_url: URL to get all links
    :return: a list of all links in the url
    """
    # create response object
    r = requests.get(archive_url, auth=(USER_NAME, PASSWORD), verify=False)
    # create beautiful-soup object
    soup = BeautifulSoup(r.content, 'html.parser')
    # find all links on web-page
    links = soup.findAll('a')
    return links


def donwload_link(link, directory):
    """

    :param link: link to donwload
    :param directory: the local directory to store it
    :return:
    """
    # obtain filename by splitting url and getting and last string
    file_name = link.split('/')[-1]

    print("Downloading file:%s" % file_name)

    # create response object
    r = requests.get(link, auth=(USER_NAME, PASSWORD), verify=False)

    # get the current working directory
    cwd = os.getcwd()
    # join the current directory with the user defined directory
    save_directory = os.path.join(cwd, directory)
    print("save directory ", save_directory)
    # creating the directory and downloading
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
    with open(save_directory + file_name, 'wb') as f:
        f.write(r.content)

    print("%s downloaded!\n" % file_name)


def download_all(archive_url, directory):
    """

    :param archive_url: The home URL
    :param directory: The home directory to store data
    :return:
    """
    urllib3.disable_warnings()
    if not os.path.exists(directory):
        os.makedirs(directory)
    links = get_links(archive_url)

    for link in links:
        if link['href'][-4:] in EXTENSION:
            donwload_link(archive_url + link['href'], directory)
        elif link['href'][0:5] in TILE:
            download_all(archive_url + link['href'], os.path.join(directory, link['href']))
        elif link['href'][0:4] in YEAR:
            download_all(archive_url + link['href'], os.path.join(directory, link['href']))
        elif link['href'][0:3] in SATELLITE:
            download_all(archive_url + link['href'], os.path.join(directory, link['href']))


def main():
    directory = os.path.join('data', 'hls_tmp')
    archive_url = "https://192.168.1.12/data-04/hls_tmp/"
    download_all(archive_url, directory)


if __name__ == "__main__":
    main()
