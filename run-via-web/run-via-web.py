from requests import Session
import json
from pathlib import Path
import wget
import sys
import argparse

DEFAULT_SETTING = {
  "is_fits_options[brightnessQuartileForASinhBeta]": "0.75",
  "is_fits_options[asinhApplications]": "2",
  "is_fits_options[brightnessQuartileForASinhAlpha]": "0.25",
  "is_fits_options[compute-starmask]": "1",
  "main_options[medFiltRad]": "1",
  "main_options[ctrDriftThresForStarMask]": "2.5",
  "main_options[useDeProjectStretch]": "1",
  "main_options[deleteClusterContainingCenter]": "1",
  "main_options[allowArcBeyond2pi]": "1",
  "main_options[errRatioThres]": "2.5",
  "main_options[mergeChkMinClusSz]": "25",
  "main_options[clusSizeCutoff]": "150",
  "main_options[stopThres]": "0.15",
  "main_options[unsharpMaskSigma]": "25",
  "main_options[unsharpMaskAmt]": "6",
  "main_options[fitUsingNonUsmIVals]": "1",
  "main_options[lookForBulge]": "1",
  "main_options[numOrientationFieldLevels]": "3",
  "main_options[useTwoStageCtrFinding]": "1",
  "main_options[useImageStandardization]": "1",
  "main_options[resizeDims]": "'[256 256]'"
}

def request_preprocess(fp):
  """
  Uploads the image and sends a request to SpArcFiRe website to start
  preprocessing the file.

  :param fp: a Path object to the source file to be uploaded and processed
  :return: response text from the request
  """
  session = Session()
  response = session.post(
    url='http://sparcfire.ics.uci.edu/process/preprocess.php',
    files={'image-selection-input': (fp.name, fp.open('rb'))},
    data=DEFAULT_SETTING
  )
  return response.text

def get_query_id(response_text):
  """
  Returns the query ID from the preprocess request.

  :param response_text: a response string from the preprocess request
  :return: a string of query ID
  """
  obj = json.loads(response_text)
  url = obj['data']['url']
  query_id = url.split('=')[1]
  return query_id

def request_process(query_id):
  """
  Sends a request to SpArcFiRe website to start processing the given query.

  :param query_id: a string of the query ID to be processed
  """
  session = Session()
  response = session.post(
    url='http://sparcfire.ics.uci.edu/process/process.php',
    data={'id': query_id},
    headers={'Content-Type': 'application/x-www-form-urlencoded; charset=UTF-8'}
  )
  assert_success(response.text)

def assert_success(response_text):
  """
  Makes sure that process request is completed successfully.

  :param response_text: a response string from the process request
  """
  obj = json.loads(response_text)
  assert obj['success'] == True

def download(destdirp, query_id, galaxy_name):
  """
  Downloads galaxy data of the given query ID as a zip file.

  :param destdirp: a Path object to the directory where the galaxy data
    will be downloaded to
  :param query_id: a string of the query ID for the galaxy
  :param galaxy_name: the name of the galaxy
  """
  url = 'http://sparcfire.ics.uci.edu/process/{}/outDir/galaxy_{}_data.zip'
  url = url.format(query_id, galaxy_name)
  wget.download(url, str(destdirp))
  print >> sys.stderr, '' # adds a newline after wget progress bar

def run_pipeline(fp, destdirp, outfp, galaxy_name):
  """
  Sends a request to SpArcFiRe website to process a file and downloads galaxy
  data as a zip file.

  :param fp: a Path object to the source file
  :param destdirp: a Path object to the directory where the galaxy data
    will be downloaded to
  :param outfp: an IO stream for writing the galaxy name and its query ID
  :param galaxy_name: the name of the galaxy
  """
  query_id = get_query_id(request_preprocess(fp))
  print query_id
  if outfp is not None:
    outfp.write(' {}\n'.format(query_id))
    outfp.flush()
  request_process(query_id)
  download(destdirp, query_id, galaxy_name)

def main(srcdir, destdir, output):
  """
  Uploads images in srcdir directory to SpArcFiRe website and downloads the
  results into destdir directory.
  Optionally, logs the query ID corresponding to each galaxy into output file.

  :param srcdir: a string of a path to source images directory
  :param destdir: a string of a path to the destination directory
  :param output: a string of a filename to write the query IDs
  """
  srcdirp = Path(srcdir)
  destdirp = Path(destdir)
  assert srcdirp.is_dir()
  assert destdirp.is_dir()

  outfp = None
  if output is not None:
    outfp = open(output, 'w')
  
  i = 1
  for fp in srcdirp.iterdir():
    if fp.is_file() and fp.suffix == '.fits':
      print '#{}: {}'.format(i, fp.stem),
      if outfp is not None:
        outfp.write(fp.stem)
      run_pipeline(fp, destdirp, outfp, fp.stem)
      i += 1

  if output is not None:
    outfp.close()

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('srcdir', help='a directory of source files')
  parser.add_argument('destdir', help='a directory for output files')
  parser.add_argument('-o', '--output',
    help='an output file for queryID corresponding to galaxy names')
  args = parser.parse_args()
  main(args.srcdir, args.destdir, args.output)
