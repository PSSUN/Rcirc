# download data from public database
# public database includes {
# plantcircbase(http://ibi.zju.edu.cn/plantcircbase/index.php),
# CIRCpedia_v2(http://www.picb.ac.cn/rnomics/circpedia/)
# }
##' @export
##' @param string one of 'ath' 'zma' 'asa' 'mm10' 'rat' 'zeb' 'fly' 'worm'
##' @param string path of downloaded file


downloadCircRNA <- function(species, out) {
  if (species == "ath") {
    download.file(
      url = "http://bis.zju.edu.cn/plantcircnet/Data/Download/ath.csv",
      destfile = out,
      method = "wget"
    )
  }
  if (species == "zma") {
    download.file(
      url = "http://bis.zju.edu.cn/plantcircnet/Data/Download/zma.csv",
      destfile = out,
      method = "wget"
    )
  }
  if (species == "osa") {
    download.file(
      url = "http://bis.zju.edu.cn/plantcircnet/Data/Download/osa.csv",
      destfile = out,
      method = "wget"
    )
  }
  if (species == "bdi") {
    download.file(
      url = "http://bis.zju.edu.cn/plantcircnet/Data/Download/bdi.csv",
      destfile = out,
      method = "wget"
    )
  }
  if (species == "gma") {
    download.file(
      url = "http://bis.zju.edu.cn/plantcircnet/Data/Download/gma.csv",
      destfile = out,
      method = "wget"
    )
  }
  if (species == "hvu") {
    download.file(
      url = "http://bis.zju.edu.cn/plantcircnet/Data/Download/hvu.csv",
      destfile = out,
      method = "wget"
    )
  }
  if (species == "tae") {
    download.file(
      url = "http://bis.zju.edu.cn/plantcircnet/Data/Download/tae.csv",
      destfile = out,
      method = "wget"
    )
  }
  if (species == "sly") {
    download.file(
      url = "http://bis.zju.edu.cn/plantcircnet/Data/Download/sly.csv",
      destfile = out,
      method = "wget"
    )
  }
  if (species == "hg38") {
    download.file(
      url = "http://www.picb.ac.cn/rnomics/circpedia/static/download_cache/human_hg38_All_circRNA.csv",
      destfile = out,
      method = "wget"
    )
  }
  if (species == "mm10") {
    download.file(
      url = "http://www.picb.ac.cn/rnomics/circpedia/static/download_cache/mouse_mm10_All_circRNA.csv",
      destfile = out,
      method = "wget"
    )
  }
  if (species == "rat") {
    download.file(
      url = "http://www.picb.ac.cn/rnomics/circpedia/static/download_cache/rat_rn6_All_circRNA.csv",
      destfile = out,
      method = "wget"
    )
  }
  if (species == "zeb") {
    download.file(
      url = "http://www.picb.ac.cn/rnomics/circpedia/static/download_cache/zebrafish_danrer10_All_circRNA.csv",
      destfile = out,
      method = "wget"
    )
  }
  if (species == "fly") {
    download.file(
      url = "http://www.picb.ac.cn/rnomics/circpedia/static/download_cache/fly_dm6_All_circRNA.csv",
      destfile = out,
      method = "wget"
    )
  }
  if (species == "worm") {
    download.file(
      url = "http://www.picb.ac.cn/rnomics/circpedia/static/download_cache/worm_ce10_All_circRNA.csv",
      destfile = out,
      method = "wget"
    )
  }
}
