about the docs
--------------

The jsonctmctree project itself is hosted on github,
and the docs are hosted at readthedocs.
I'm not using the full integration option,
so the project accounts on github and readthedocs are not as tightly
linked as they would like.
This means that the github project doesn't have readthedocs flair,
and that the readthedocs docs are not automatically updated using
commit magic when changes are pushed to the github project.

The docs are in a specific docs top-level directory in the github
project, which I think is required by readthedocs.
I'm using a manually set admin option to tell readthedocs that I'm using
sphinx, otherwise readthedocs tries to use MkDocs.

I'm currently trying to write some sphinx docs following
https://docs.readthedocs.org/en/latest/getting_started.html
with some help from
http://matplotlib.org/sampledoc/getting_started.html
which unfortunately seems a little bit outdated.
As you can tell from these urls, I'm indeed just "getting started."

I am trying to follow
http://stackoverflow.com/questions/4259105
to make the json highlighting work.
Sphinx will already correctly highlight literal inclusions of
entire .json files, but it is not highlighting code-block json.
*Update*:
So it turns out that I was just giving it an invalid json
block that looked kind of like ``{1:2}`` whereas only strings
are allowed as keys.
