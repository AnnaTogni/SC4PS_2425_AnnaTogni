# HomeWork 01: Setting up a Linux machine on CloudVeneto for C development

This guide goes from zero to a working Linux virtual machine (VM) on
[CloudVeneto](https://cloudveneto.ict.unipd.it/dashboard), ready to edit, compile and
run C programs.

The VM runs **AlmaLinux 9**, the distribution used in the course. It is a free
rebuild of Red Hat Enterprise Linux and is also used at CERN. Its package
manager is `dnf` (`yum` still works as an alias).

---

## Contents
1. [Prerequisites](#1-prerequisites)
2. [Create an SSH key pair](#2-create-an-ssh-key-pair)
3. [Launch the virtual machine](#3-launch-the-virtual-machine)
4. [Connect via SSH](#4-connect-via-ssh)
5. [First steps on the VM: users, updates, packages](#5-first-steps-on-the-vm)
6. [Install the C development tools](#6-install-the-c-development-tools)
7. [Write, compile and run a C program](#7-write-compile-and-run-a-c-program)
8. [Automate the build with a Makefile](#8-automate-the-build-with-a-makefile)
9. [Moving files and using Git](#9-moving-files-and-using-git)
10. [Linux cheat-sheet](#10-linux-cheat-sheet)
11. [Troubleshooting](#11-troubleshooting)

---

## 1. Prerequisites

| What | Why |
|---|---|
| A CloudVeneto account, enrolled in the course project | To create VMs on the dashboard |
| A terminal with `ssh` | Linux/macOS: the built-in terminal. Windows: PowerShell or WSL. |
| A University network connection (or VPN), or the CloudVeneto gate | VMs have **private** IP addresses (see §4) |

Throughout the guide:
- Lines starting with `$` are typed on **your own computer**.
- Lines starting with `[vm]$` are typed **on the VM**.
- Don't type the `$` itself.

---

## 2. Create an SSH key pair

SSH logs in with a **key pair** instead of a password:
- the **private key** stays on your computer and is never shared;
- the **public key** is copied onto the VM.

There are two ways to get one.

**Option A: let CloudVeneto generate it** (the way shown in class)
1. On the dashboard go to **Project → Compute → Key Pairs → Create Key Pair**.
2. Give it a name (e.g. `sc4ps-key`) and choose type **SSH Key**.
3. The browser downloads `sc4ps-key.pem`. That is your private key, and it **cannot be downloaded again**.
4. Move it to `~/.ssh` and protect it. SSH refuses keys that other users can read:
   ```bash
   $ mv ~/Downloads/sc4ps-key.pem ~/.ssh/
   $ chmod 600 ~/.ssh/sc4ps-key.pem
   ```

**Option B: generate it yourself and import the public part**
```bash
$ ssh-keygen -t ed25519 -f ~/.ssh/sc4ps-key   # creates sc4ps-key (private) + sc4ps-key.pub (public)
$ cat ~/.ssh/sc4ps-key.pub                     # copy this line
```
Then on the dashboard click **Key Pairs → Import Public Key** and paste the line.

---

## 3. Launch the virtual machine

1. **Allow SSH traffic.** Go to **Network → Security Groups**, open the group you will use (or create
   `ssh-only`), click **Add Rule**, select **SSH** (TCP port 22) and confirm.
2. **Create the instance.** Go to **Compute → Instances → Launch Instance** and fill in the tabs:

   | Tab | Setting |
   |---|---|
   | Details | Instance name, e.g. `sc4ps-vm` |
   | Source | Boot source **Image**, pick an **AlmaLinux 9** image |
   | Flavor | A small one is enough for this course (e.g. 2 vCPU / 4 GB RAM) |
   | Networks | The project's network (usually pre-selected) |
   | Security Groups | `default` + the group with the SSH rule |
   | Key Pair | `sc4ps-key` |

3. Click **Launch Instance**. Wait until the status becomes **Active**, then write down the
   VM's **IP address** (e.g. `10.67.x.x`).

> Each image's description on the dashboard states its **default login user**. For
> AlmaLinux images this is typically `almalinux`. Use that name below in place of `<user>`.

---

## 4. Connect via SSH

**From inside the University network (or VPN)** connect directly:
```bash
$ ssh -i ~/.ssh/sc4ps-key.pem <user>@<vm-ip>
```

**From outside**, the VM's private IP is not reachable. You have to "jump" through the CloudVeneto
gate machine, using your CloudVeneto username:
```bash
$ ssh -J <cloudveneto-username>@gate.cloudveneto.it -i ~/.ssh/sc4ps-key.pem <user>@<vm-ip>
```

The first time, SSH asks *"Are you sure you want to continue connecting?"*. Type `yes`.
If you see the VM's prompt, e.g. `[almalinux@sc4ps-vm ~]$`, you're in!

### Make it a one-word command
Add this to `~/.ssh/config` on your computer:
```
Host cvgate
    HostName gate.cloudveneto.it
    User <cloudveneto-username>

Host sc4ps
    HostName <vm-ip>
    User <user>
    IdentityFile ~/.ssh/sc4ps-key.pem
    ProxyJump cvgate        # remove this line when you are on the University network
```
From now on, `ssh sc4ps` is enough.

> **Tip:** VS Code's **Remote - SSH** extension uses the same `~/.ssh/config`. Open the command
> palette, choose *Remote-SSH: Connect to Host…* and select `sc4ps`. You can then edit files on the
> VM with a full editor.

---

## 5. First steps on the VM

```bash
[vm]$ cat /etc/os-release        # check: AlmaLinux 9
[vm]$ sudo dnf -y update         # update all installed packages
```

**`sudo`** ("super-user do") runs a single command with administrator rights. You need it to
install software. If it is missing, install it as root with `dnf install sudo`.

Useful `dnf` commands:

| Command | Effect |
|---|---|
| `dnf search <word>` | Search for a package |
| `sudo dnf install <pkg>` | Install a package |
| `sudo dnf remove <pkg>` | Uninstall it |
| `dnf list installed` | Show what is installed |
| `sudo dnf update` | Update everything |

---

## 6. Install the C development tools

The **"Development Tools"** group installs everything needed to build C code: the `gcc` compiler,
`make`, `gdb` (debugger), `git` and more.

```bash
[vm]$ sudo dnf -y groupinstall "Development Tools"
```

Some useful extra packages live in the **CRB** (CodeReady Builder) and **EPEL** repositories.
Enable them once:

```bash
[vm]$ sudo dnf -y install dnf-plugins-core epel-release
[vm]$ sudo dnf config-manager --set-enabled crb
```

Optional but handy:

```bash
[vm]$ sudo dnf -y install vim emacs-nox htop valgrind
```

Check that everything is there:

```bash
[vm]$ gcc --version
[vm]$ make --version
[vm]$ gdb --version
```

---

## 7. Write, compile and run a C program

### 7.1 Create a project folder and a source file
```bash
[vm]$ mkdir -p ~/sc4ps/hello && cd ~/sc4ps/hello
[vm]$ vim hello.c
```

In **Vim**:
- Press `i` to start typing (*insert mode*) and write the program below.
- Press `Esc`, type `:wq` and press Enter to **w**rite and **q**uit.

In **Emacs** (`emacs -nw hello.c`):
- Just type.
- `Ctrl+X Ctrl+S` saves; `Ctrl+X Ctrl+C` quits.

```c
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("Hello from CloudVeneto!\n");
    printf("sqrt(2) = %.15f\n", sqrt(2.0));
    return 0;
}
```

### 7.2 Compile
```bash
[vm]$ gcc -Wall -Wextra -O2 -o hello hello.c -lm
```

| Piece | Meaning |
|---|---|
| `gcc` | The GNU C compiler |
| `-Wall -Wextra` | Turn on warnings. Always read them, they catch real bugs. |
| `-O2` | Optimise the code (use `-O0 -g` when debugging with `gdb`) |
| `-o hello` | Name of the output executable (default would be `a.out`) |
| `hello.c` | The source file |
| `-lm` | Link the **m**ath library (needed for `sqrt`, `sin`, `exp`, …) |

### 7.3 Run
```bash
[vm]$ ./hello
Hello from CloudVeneto!
sqrt(2) = 1.414213562373095
```

The `./` means "the program in *this* directory". Linux doesn't look in the current folder by
default.

### 7.4 Debugging in short
```bash
[vm]$ gcc -Wall -O0 -g -o hello hello.c -lm   # -g keeps debug information
[vm]$ gdb ./hello                             # break main, run, next, print x, quit
[vm]$ valgrind ./hello                        # detects memory errors and leaks
```

---

## 8. Automate the build with a Makefile

Typing the `gcc` line every time is error-prone. A `Makefile` stores the build recipe, and
`make` rebuilds only what changed.

```makefile
CC     = gcc
CFLAGS = -Wall -Wextra -O2
LDLIBS = -lm

hello: hello.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f hello

.PHONY: clean
```

> The indented lines **must start with a TAB**, not spaces.

```bash
[vm]$ make          # builds ./hello
[vm]$ make clean    # removes it
```

---

## 9. Moving files and using Git

**Copy files between your computer and the VM** with `scp`. It uses the same `~/.ssh/config`
alias as `ssh`:
```bash
$ scp hello.c sc4ps:~/sc4ps/hello/      # computer → VM
$ scp sc4ps:~/sc4ps/hello/out.txt .     # VM → computer
```

**Git** (installed with the Development Tools) keeps the homework under version control.
To use it, configure your identity once and clone your repository:
```bash
[vm]$ git config --global user.name  "Your Name"
[vm]$ git config --global user.email "you@example.com"
[vm]$ git clone git@github.com:<github-user>/<repo>.git
```

Cloning over SSH requires a key on the VM:
1. Create one with `ssh-keygen -t ed25519`.
2. Add `~/.ssh/id_ed25519.pub` on GitHub under **Settings → SSH keys**.

Then the everyday cycle is `git add .` → `git commit -m "message"` → `git push`.

---

## 10. Linux cheat-sheet

| Command | What it does |
|---|---|
| `pwd` | Print the current directory |
| `ls -lah` | List all files, with sizes and permissions |
| `cd <dir>` / `cd ..` / `cd ~` | Change directory / go up one / go home |
| `mkdir -p a/b` | Create directories |
| `cp src dst`, `mv old new`, `rm file` | Copy, move/rename, delete (**no recycle bin!**) |
| `cat`, `less`, `head`, `tail` | Show file contents |
| `grep "word" file` | Search inside files |
| `chmod +x script.sh` | Make a file executable |
| `cat /proc/cpuinfo`, `nproc` | CPU information, number of cores |
| `free -h`, `cat /proc/meminfo` | Memory (RAM) usage |
| `df -h`, `du -sh *` | Disk space, size of each item here |
| `top` / `htop` | Running processes (`q` to quit) |
| `man <cmd>` / `<cmd> --help` | Documentation |
| `exit` or `Ctrl+D` | Log out of the VM |

**Small scripting exercise** to put it all together:
```bash
[vm]$ touch myscript.sh
[vm]$ vim myscript.sh        # write: echo "I have $(nproc) cores"
[vm]$ chmod +x myscript.sh
[vm]$ ./myscript.sh
```

---

## 11. Troubleshooting

| Problem | Fix |
|---|---|
| `WARNING: UNPROTECTED PRIVATE KEY FILE!` | `chmod 600 ~/.ssh/sc4ps-key.pem` |
| `Permission denied (publickey)` | Wrong user name or key. Check the image's default user and the `-i` path. |
| `Connection timed out` | You are outside the University network: use the gate (`-J`). Check the security group allows port 22. |
| `gcc: command not found` | Development Tools not installed (§6) |
| `undefined reference to 'sqrt'` | Add `-lm` at the **end** of the `gcc` line |
| `make: *** missing separator` | The Makefile recipe lines start with spaces instead of a TAB |
| `bash: ./hello: Permission denied` | `chmod +x hello`, or recompile |
| Disk full | `df -h` to check, `du -sh *` to find what is big |

---

### References
- CloudVeneto user guide: <https://userguide.cloudveneto.it>
- AlmaLinux documentation: <https://wiki.almalinux.org>
- Linux tutorials: <https://github.com/cosmicqbit/Learn-Linux>, <https://github.com/labex-labs/linux-free-tutorials>
- GCC options: `man gcc`, or <https://gcc.gnu.org/onlinedocs/>
