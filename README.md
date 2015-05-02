# Cluster Install / Management Scripts

A few scripts I use to setup and maintain an [Arch
Linux](https://www.archlinux.org/) compute cluster.
Probably some bugs that need to be worked out, as I keep adding to the install
scripts as I modify our newest cluster. I got tired of having to install a newer
compiler or a newer version of Python every time a new version of software we
needed came out. Arch Linux is a rolling release so that issue is solved. Also I
only have to install what I need (I don't really need Gnome for this). These
scripts help me to keep track what I've done, as well as reproduce the
environment in the future.

## Install Scripts

Install scripts are for use with an Arch Linux live image. Download the image
and boot to the USB drive as normal. The scripts follow the [Beginner's
Guide](https://wiki.archlinux.org/index.php/Beginners'_guide)
plus some additional installations and modifications. The scripts are used both
on the master and other nodes (user is prompted to select which one when
appropriate). When in the live environment simply
download the first script:

    wget https://raw.githubusercontent.com/wesbarnett/cluster/master

Or with a URL shortener:

    wget http://bit.ly/tulane-install

### install

First make the script executable. If you used the URL shortener then the script
will be named `tulane-install` instead of just `install`. Then run the script in
the live environment.

This script:

* Partitions the disk(s) (currently setup for two disks)
* Creates the file system
* Mounts the partitions
* Installs the base system
* Generates fstab
* Grabs part 2, makes it executable, and puts it in /mnt
* Switches to chroot environment

### install-part2

After you are in the chroot environment you should see the second script. Run
the script.

This script:

* Sets the locale (English)
* Sets the time zone (Central US)
* Sets the hostname ("master" for the master node, "nodeX" where X is node
  number for the nodes)
* Sets up networking with
  [systemd-networkd](https://wiki.archlinux.org/index.php/Systemd-networkd) (Master is static ip 192.168.0.1,
  nodes are all 192.168.0.X, and IP forwarding is set up on the master).
* Updates /etc/hosts with nodes
* Creates initial ramdisk
* Sets root password
* Install [GRUB](https://wiki.archlinux.org/index.php/GRUB) as bootloader
* Sets up [NFS](https://wiki.archlinux.org/index.php/NFS). Specifically the
  [pacman cache on the master will be
shared](https://wiki.archlinux.org/index.php/Pacman_tips#Network_shared_pacman_cache)
with the nodes. Also the home directories on the nodes will be shared with the
master.
* Downloads part 3, and makes executable.

### install-part3

After part2 is down, you plug the ethernet cable of the node (if not on master)
into the internal switch now, since the packages can be retrieved from the
master's cache. After you've restarted you should see part 3 already downloaded.

This script:

* Downloads and installs: `vim`, `openssh`, `cmake`, `git`, `mlocate`, `zsh`, `screen`,
  `xorg-server`, `fftw`, `xorg-auth`
* Updates the locate database.
* [Enables the `wheel` group in
  `sudoers`](https://wiki.archlinux.org/index.php/Sudo).
* Adds your user (be sure to change the variable at the beginning of the
  script).
* [Adds all other
  users](https://wiki.archlinux.org/index.php/Users_and_groups#Example_adding_a_user),
setting a default password for them, and adds a /data/USER directory for the
second disk (assumes you have two disks). Be sure
to change the default password and user names at the beginning of the script.
* Sets up iptables [simple stateful
  firewall](https://wiki.archlinux.org/index.php/Simple_stateful_firewall) for master, enables [sharing of the
  internet connection](https://wiki.archlinux.org/index.php/Internet_sharing) with the nodes, and enables [sshguard](https://wiki.archlinux.org/index.php/Sshguard) on the master to stop ssh attacks.
* Sets up [sshd](https://wiki.archlinux.org/index.php/Secure_Shell), disabling
  root login and enabling X forwarding.
* Sets up [pdnsd](https://wiki.archlinux.org/index.php/Pdnsd) for DNS cache
  sharing on the master with the nodes.
* [Syncs time with
  NTP](https://wiki.archlinux.org/index.php/Time#Time_synchronization).
* Installs [Intel microcode
  update](https://wiki.archlinux.org/index.php/Microcode).
* Enables
  [auto-logout](https://wiki.archlinux.org/index.php/Security#Automatic_logout) for virtual consoles.

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
