import ndex.client as nc
import io
import requests

source_names = ["CALM3"]
target_names = ["NFATC2"]
host = "http://www.ndexbio.org"
directed_path_query_url = "http://general.bigmech.ndexbio.org:5603/directedpath/query"

ndex = nc.Ndex(host=host)

reference_network_cx = ndex.get_network_as_cx_stream("5294f70b-618f-11e5-8ac5-06603eb7f303")

#====================
# Assemble REST url
#====================
target = ",".join(target_names)
source = ",".join(source_names)
max_number_of_paths = 20
url = directed_path_query_url + '?source=' + source + '&target=' + target + '&pathnum=' + str(max_number_of_paths)

f = io.BytesIO()
f.write(reference_network_cx.content)
f.seek(0)

url = directed_path_query_url + '?source=' + source + '&target=' + target + '&uuid=5294f70b-618f-11e5-8ac5-06603eb7f303&server=www.ndexbio.org&pathnum=' + str(max_number_of_paths)

r = requests.post(url)

print r.content
