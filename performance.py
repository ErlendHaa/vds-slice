import click
import time
import json
import requests

import numpy as np

from pprint import pformat
from requests_toolbelt.multipart.decoder import MultipartDecoder

# Highly assumes that you run the servers on the respective ports
hosts = {
    'master'              : 'http://localhost:8080',
    'go-single-handle'    : 'http://localhost:8081',
    'go-multiple-handle'  : 'http://localhost:8082',
    'cpp-single-handle'   : 'http://localhost:8083',
    'cpp-multiple-handle' : 'http://localhost:8084',
    'grpc'                : 'http://localhost:8090',
    'radix-dev'           : 'https://server-oneseismictest-dev.playground.radix.equinor.com',
    'radix-test'          : 'https://server-oneseismictest-test.playground.radix.equinor.com',
    'radix-grpc'          : 'https://scheduler-oneseismictest-grpc.playground.radix.equinor.com'
}

class OneseismicClient:
    def __init__(self, host, vds, sas = None):
        self.host = host
        self.vds  = vds
        self.sas  = sas

    def metadata(self):
        """Get metadata from the OpenVDS volume"""

        response = requests.post(f'{self.host}/metadata',
            headers = { 'Content-Type' : 'application/json' },
            data    = json.dumps({
                'vds' : self.vds,
                'sas' : self.sas
            })
        )
        if not response.ok: raise RuntimeError(response.text)
        return response.json()

    def fence(self, coordinates, system):
        response = requests.post(f'{self.host}/fence',
            headers = {'Content-Type': 'application/json'},
            data = json.dumps({
                'coordinates'      : coordinates,
                'coordinateSystem' : system,
                'vds'              : self.vds,
                'sas'              : self.sas
            })
        )
        if not response.ok: raise RuntimeError(response.text)
        parts = MultipartDecoder.from_response(response).parts

        meta = json.loads(parts[0].content)
        data = parts[1].content

        return np.ndarray(meta['shape'], meta['format'], data), meta

def send(host, vds, sas, *args, log=True, runs = 5, equalto = None):
    hostname = list(hosts.keys())[list(hosts.values()).index(host)]
    if log:
        print(f'Fetching fence on host: "{hostname}" ({host})')


    client = OneseismicClient(host, vds, sas)

    durations = []
    for _ in range(runs):
        start = time.perf_counter()
        data, meta = client.fence(*args)
        stop = time.perf_counter()

        durations.append(stop - start)

    if equalto is not None:
        np.testing.assert_array_equal(data, equalto)

    if log:
        size = data.nbytes / (1024 * 1024)
        average = sum(durations) / runs
        minimum = min(durations)
        maximum = max(durations)

        print(f'Fetching fence {runs} time(s). Fence size: {size:.3f} MB')
        print(f'Avg   | Min   | Max')
        print(f'------|-------|------')
        print(f'{average:.3f} | {minimum:.3f} | {maximum:.3f}\n')

    return data, meta

description = """
This program expect the following implementations runs on these endpoints
(sorry for the formatting. I blame click):

{}""".format(json.dumps(hosts, indent = 4))

@click.command(help = description)
@click.option('--sas', required = True, type = str, help = "sas-token for storage account 'vdsbenchmark'")
def main(sas):
    vds = 'https://vdsbenchmark.blob.core.windows.net/gigamerge/100_Merge_surveys_S1_S10/A6C1BD8CC0C7F460'
    
    # This number means that we can split this request in 2,4,8,16,32 parts
    # without having multiple parts having to download the same chunk. I.e. the
    # list of coordinates is always split on chunk boarder.
    # This is deliberate to counter measure to our stupid splitting logic that
    # would otherwise result in double work.
    coordinates = [(x, 5064) for x in range(16384)]

    ## Warm up blob store
    d0, _ = send(hosts['master'], vds, sas, coordinates, 'ij', log = False, runs = 1)

    send(hosts['master'],              vds, sas, coordinates, 'ij', equalto = d0)
    send(hosts['go-single-handle'],    vds, sas, coordinates, 'ij', equalto = d0)
    send(hosts['go-multiple-handle'],  vds, sas, coordinates, 'ij', equalto = d0)
    send(hosts['cpp-single-handle'],   vds, sas, coordinates, 'ij', equalto = d0)
    send(hosts['cpp-multiple-handle'], vds, sas, coordinates, 'ij', equalto = d0)
    send(hosts['grpc'],                vds, sas, coordinates, 'ij', equalto = d0)
    ## TODO: radix-dev and radix-test has to little memory to serve this request
    # send(hosts['radix-dev'],           vds, sas, coordinates, 'ij', equalto = d0)
    # send(hosts['radix-test'],          vds, sas, coordinates, 'ij', equalto = d0)
    send(hosts['radix-grpc'],          vds, sas, coordinates, 'ij', equalto = d0)


if __name__ == '__main__':
    main()
