# hIPPYlib / web

This repo contains the [MkDocs](http://mkdocs.org) source files of the hIPPYlib website.

## Prerequisites

* Install [MkDocs](http://mkdocs.org)
    sudo pip install mkdocs==0.14.0

* Install [bootstrap](http://getbootstrap.com/) and [bootswatch](https://bootswatch.com/) themes
    sudo pip install mkdocs-bootstrap==0.1.1
    sudo pip install mkdocs-bootswatch==0.1.0
    
* Install [MathJax](https://www.mathjax.org/) support using [python-markdown-math](https://github.com/mitya57/python-markdown-math)
    sudo pip install python-markdown-math==0.2


## To make changes to the website

* Clone this repository
    git clone git@bitbucket.org:hippylibdev/hippylib-website.git
* Add the github repository as public (i.e. where we will publish the website).
    git remote add public https://github.com/hippylib/hippylib.github.io.git
* Edit or add some `.md` files in the src folder (you may also need to update the `mkdocs.yml` config); 
* Preview locally
    mkdocs serve
* Publish on GitHub
    git push origin master:master
    mkdocs gh-deploy --clean --remote-branch master --remote-name public
    git reset --hard origin/master
