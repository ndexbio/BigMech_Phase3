__author__ = 'danielcarlin'

import argparse
import ndex.client as nc
import ndex.networkn as networkn

net_uuid='b04e406b-dc88-11e6-86b1-0ac135e8bacf'

G = networkn.NdexGraph(server='http://public.ndexbio.org', uuid=net_uuid)
