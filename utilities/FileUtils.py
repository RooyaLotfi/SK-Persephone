import re


class FileUtils:
    def __init__(self):
        pass

    def atoi(self, text):
        return int(text) if text.isdigit() else text

    def natural_keys(self, text):
        return [self.atoi(c) for c in re.split(r'(\d+)', text)]

    def parse_tile_file(self, image):
        file_name = image.replace('.', ' ').replace('_', ' ').split()

        if file_name[0] == 'HLS':
            sat_id = file_name[1]
            tile = file_name[2]
            time = file_name[3]

            for i in range(0, len(time), 4):
                year = time[:4]
                day = time[i:i + 4]
            band_no = file_name[8]

        elif file_name[0] == 'MCD43A4':
            sat_id = file_name[0]
            year = file_name[1][1:5]
            day = file_name[1][5:]
            tile = file_name[2]
            band_no = file_name[5]

        return sat_id, tile, year, day, band_no

    def parse_mask_file(self, image):
        file_name = image.replace('.', '_').split('_')
        tile = file_name[2]
        year = file_name[3]
        crop_type = file_name[4]

        return tile, year, crop_type

    def parse_qamasked_file(self, image):
        file_name = image.replace('.', '_').split('_')
        sat_id = file_name[0]
        tile = file_name[1]
        year = file_name[2]
        day = file_name[3]
        band_no = file_name[4]

        return sat_id, tile, year, day, band_no

    def parse_preprocessed_file(self, image):
        file_name = image.replace('.', '_').split('_')
        sat_id = file_name[0]
        tile = file_name[1]
        year = file_name[2]
        day = file_name[3]
        band_no = file_name[4]
        crop_type = file_name[5]

        return sat_id, tile, year, day, band_no, crop_type

    def parse_vi_file(self, text_file):
        file_name = text_file.replace('.', '_').split('_')
        tile = file_name[0]
        sat_id = file_name[1]
        year = file_name[2]
        day = file_name[3]
        crop_type = file_name[4]
        vi = file_name[5]

        return sat_id, tile, year, day, crop_type, vi

    def parse_landcover_file(self, path):
        path = path.replace('.', '_').replace('/', '_').split('_')
        tile = path[4]
        sat_id = path[5]
        year = path[6]

        return sat_id, tile, year

    def get_days(self, file_lists):
        days = []
        for files in file_lists:
            for file in files:
                if ".txt" in file or ".tif" in file:
                    file_name = file.replace('.', '_').split('_')
                    day = file_name[3]
                    days.append(day)
        days = list(set(days))

        return days

    def get_years(self, file_lists):
        years = []
        for files in file_lists:
            for file in files:
                if ".txt" in file or ".tif" in file:
                    file_name = file.replace('.', '_').split('_')
                    year = file_name[2]
                    years.append(year)
        years = list(set(years))

        return years

    def get_tiles(self, file_lists):
        tiles = []
        for files in file_lists:
            for file in files:
                if ".txt" in file or ".tif" in file:
                    file_name = file.replace('.', '_').split('_')
                    tile = file_name[1]
                    tiles.append(tile)
        tiles = list(set(tiles))

        return tiles

    def parse_csv_file(self, csv_file):
        file_name = csv_file.replace('.', '_').split('_')
        day = file_name[1]

        return day

    def parse_data_file(self, data_file):
        file_name = data_file.replace('.', '_').split('_')
        year = file_name[1]
        day = file_name[2]

        return year, day