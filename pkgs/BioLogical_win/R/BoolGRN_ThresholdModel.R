#' @title 8 threshold-based genetic networks
#' @description A list containing networks from the papers, with each element 
#' representing an independent network. Their names briefly describe the features 
#' of the networks. 
#' @details The source of the network is provided in the references:
#' \itemize{
#'   \item \code{ArabidopsisFlower}: The network determines during Arabidopsis thaliana flower development \href{https://doi.org/10.1105/tpc.104.021725}{Espinosa-Soto2004}
#'   \item \code{emt26network}: 36-node epithelial-to-mesenchymal transition (EMT) network \href{https://doi.org/10.1038/npjsba.2015.14}{Steinway2015}
#'   \item \code{emt72network}: 72-node EMT network \href{https://doi.org/10.1088/1478-3975/aaf8d4}{Jia2019}
#'   \item \code{gastricnetwork}: Endogenous molecular network of gastric cancer \href{https://doi.org/10.18632/oncotarget.3633}{Li2015}
#'   \item \code{sclcnetwork}: Transcription factor network in small cell lung cancer \href{https://doi.org/10.1158/0008-5472.CAN-16-1467}{Udyavar2016}
#'   \item \code{stemcellnetwork}: Network of induced pluripotent stem cells \href{https://doi.org/10.1371/journal.pcbi.1002300}{Chang2011}
#'   \item \code{yeastbudding}: Budding yeast cell cycle network \href{https://doi.org/10.1073/pnas.0305937101}{Li2004}
#'   \item \code{yeastfission}: Fission yeast cell cycle network \href{https://doi.org/10.1371/journal.pone.0001672}{Davidich2008}
#' }
#' Each element is a four-element formatted list that records information of 
#' a Boolean network (See \link{BoolGRN_CellCollective}). Other detail information 
#' about genetic networks, please see original papers.
#' @usage
#' data(BoolGRN_ThresholdModel)
"BoolGRN_ThresholdModel"
