'''
Modified by mor.chalabi
This API is provided by Enrichr: https://maayanlab.cloud/Enrichr/help#api
'''

import json
import requests

def Enrichr(gene_set, lib_):

  # querying
  
  ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'     # URL to send our gene list to
  
  with open(gene_set, 'r') as file:                           # reading in gene set
      genes_str = file.read()
  
  payload = { 'list': (None, genes_str) }
  
  response = requests.post(ENRICHR_URL, files = payload)      # sending our request
  if not response.ok:
      raise Exception('Error analyzing gene list')
  
  resp_ = json.loads(response.text)                           # result of request; it contains our request's userListId created by Enrichr web server
  
  # enrichment analysis
  
  ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
  query_string = '?userListId=%s&backgroundType=%s'
  user_list_id = resp_.get('userListId')
  gene_set_library = lib_                                     # library to query against
  
  response = requests.get(
      ENRICHR_URL + query_string % (user_list_id, gene_set_library)
   )
  if not response.ok:
      raise Exception('Error fetching enrichment results')
  
  return json.loads(response.text)
