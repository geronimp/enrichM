(define-module (ace packages enrichm)
  #:use-module (gnu packages)
  #:use-module (gnu packages bioinformatics)
  #:use-module (gnu packages python)
  #:use-module (gnu packages statistics)
  #:use-module (gnu packages time)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix packages)
  #:use-module (guix download)
  #:use-module (guix git-download)
  #:use-module (guix gexp)
  #:use-module (guix utils)
  #:use-module (guix build-system gnu)
  #:use-module (guix build-system python)
  #:use-module (gnu packages algebra)
  #:use-module (gnu packages bioinformatics)
  #:use-module (gnu packages boost)
  #:use-module (gnu packages check)
  #:use-module (gnu packages compression)
  #:use-module (gnu packages databases)
  #:use-module (gnu packages python)
  #:use-module (gnu packages xml)
  #:use-module (gnu packages machine-learning)
  #:use-module (gnu packages web))

(define-public enrichm
  (let ((commit "8fbb39f7ce9c95aba9fbcff6ac157dd9e1a0c2f4"))
    (package
     (name "enrichm")
     (version "0.4.6")
     (source (origin
              (method git-fetch)
              (uri (git-reference
                    (url "https://github.com/geronimp/enrichM.git")
                    (commit commit)))
              (file-name (string-append name "-" version "-checkout"))
              (sha256
               (base32
                "1y83l5ywv124292vk5agpyjlhr0l5bl6qf9czjg525bi5m6imy1k"))))
    (build-system python-build-system)
    (arguments
     `(#:python ,python ; python-3 only
       ))
    (inputs
     `(("python-dateutil" ,python-dateutil)
       ("python-statsmodels" ,python-statsmodels)
       ("python-numpy" ,python-numpy)
       ("python-pandas" ,python-pandas)
       ("python-scipy" ,python-scipy)
       ("python-biopython" ,python-biopython)
       ("python-six" ,python-six)
       ("python-scikit-learn" ,python-scikit-learn)))
    (home-page "https://github.com/geronimp/enrichM")
    (synopsis "A toolbox to compare functional composition of population genomes")
    (description "EnrichM is a toolbox for comparing the functional composition of population genomes.")
    (license license:gpl3+))))
enrichm