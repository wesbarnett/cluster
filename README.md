# Cluster Management Scripts

A few scripts I use to maintain an [Arch
Linux](https://www.archlinux.org/) compute cluster.
Probably some bugs that need to be worked out, as I keep adding to the install
scripts as I modify our newest cluster. I got tired of having to install a newer
compiler or a newer version of Python every time a new version of software we
needed came out. Arch Linux is a rolling release so that issue is solved. Also I
only have to install what I need (I don't really need Gnome for this). These
scripts help me to keep track what I've done, as well as reproduce the
environment in the future.

## Pkgbuilds

These are some [PKGBUILD](https://wiki.archlinux.org/index.php/PKGBUILD) scripts
for installing some software not found in the repository. The package only needs
to be build once which can then be copied to
all of the nodes and installed. This system allows
[pacman](https://wiki.archlinux.org/index.php/Pacman) to track what has been
installed, so packages can easily be upgraded or uninstalled.

## Management Scripts

Currently consists of two [tmux](https://wiki.archlinux.org/index.php/Tmux)
scripts which allow one to open up a terminal with a window to every node and
type one command into allow nodes at the same
time. Useful for installing / updating all nodes at the same time.
