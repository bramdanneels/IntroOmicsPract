# Shell tips & tricks

This document can act as a reference of some often used shell commands that can make life easier when working on the command line.

## Navigating the filesystem

| Command | Usage |
| --- | --- |
| `cd` | Change directory; Go to home directory |
| `cd ..` | Go to previous (parent) directory |
| `cd /specific/directory` | Go to the `/specific/directory` |
| `pwd` | Print current working directory to the screen |
| `ls` | List files in current directory |
| `ls -lah` | List all files, including hidden (`-a`), in long format (`-l`), using human-readable sizes (`-h`) |
| `ls /some/folder` | List the content of `/some/folder` |

## Reading files

| Command | Usage |
| --- | --- |
| `cat file1.txt` | Print the contents of `file1.txt` to the screen |
| `head -n 50 file2.txt` | Print the first 50 lines of `file2.txt` to the screen |
| `tail -n 50 file3.txt` | Print the last 50 lines of `file3.txt` to the screen |
| `less file4.txt` | Read the content of `file4.txt` in the `less` program |

## Manipulating files and folders

| Command | Usage |
| --- | --- |
| `mkdir newfolder` | Create a new folder in the current directory, called `newfolder` |
| `mkdir -p /create/new` | Create a folder called `new` in the folder `/create`, and create the parent folder(s) (`/create`) if they don't exist. `-p` will also not give an error if the folder you are trying to create already exists |
| `touch file.txt` | Create an empty file named `file.txt` |
| `cat > file.txt` | Start writing to file. Everything typed into the terminal will be stored in `file.txt`. Press `Ctr+C` to stop writing |
| `cat "some command" > file.txt` | Write the string "some command" to `file.txt` |
| `cat fileA.txt fileB.txt > fileAB.txt` | Write the contents of fileA and fileB to fileAB |
| `mv fileA.txt folder/fileA.txt` | Move fileA to folder/ |
| `mv fiel1.txt fileA.txt` | Rename fielA.txt to fileA.txt |
| `cp file1.txt file2.txt` | Copy file1 to file2 |
| `ln -s /full/path/to/file1 file2` | Make a symbolic link from file1 to file2. This means that if a program wants to read file2, it will be redirected to file1 isntead. If file2 is modified, a copy with the modifications will be made |
| `ln /full/path/to/file1 file2` | Make a hard link from file1 to file2. This means that if a program wants to read file2, it will be redirected to file1 instead. If file2 is modified, file1 is modified instead |
| `nano file1.txt` | Allows modifying file1 in the command line text editor `nano`. Other editors exits as well (`vim`,  |

## Miscellaneous commands

### Grep, awk, & sed

Grep, awk, and sed are three of the most usefull programs in linux. They allow you to search files (`grep`), make calculations (`awk`), and replace patterns (`sed`).
Below you can find some examples of each, but this is just the tip of the iceberg regarding their functionalities.

#### Grep

Grep is used to search for the presence of patterns in a text.

| Command | Usage |
| --- | --- |
| `grep "pattern" file.txt` | Program to search for patterns in text. This command will search and print all lines wich contain "pattern" in the file "file.txt" |
| `grep -c "pattern" file.txt` | Counts in how many lines "pattern" is present in file.txt |
| `grep -v "pattern" file.txt` | Print all lines which do not contain "pattern"  |
| `grep -f patterns.lst file.txt` | Print the lines which contain any pattern of the patterns stored in "patterns.lst" |
| `grep -i "pattern" file.txt` | Print lines containing "pattern, but ignore the case (e.g. "PaTTerN" will count as a match) |
| `grep -w "pattern" file.txt` | Only print lines that have the word "pattern" (e.g. "patterns" will not count as a match) |

> Many more options exist, see `man grep` for a full overview.
> Grep can also use regular expressions (specific patterns that match) using the `-e` option. 
> For example: `grep -e "^Str[io]ng"` will print lines which start with either "Strong" or "String".

Sed is used to replace patterns in a text.
Sed generally works as `sed 's/pattern/replacement/'`, where it will replace "pattern" by "replacement".
Standardly, `sed` will only replace the first instance, add a `g` after the last slash to do it over the whole file (globally).

| Command | Usage |
| --- | --- |
| `sed 's/human/monkey/g` file.txt | Print the content of "file.txt", but replace all instances of "human" in a tekst by "monkey" |
| `sed 's/A/a/g' file.txt` | Print the content of "file.txt", but replace all uppercase As with lowercase as |
| `sed 's/first/last/' file.txt` | Print the content of "file.txt", but replace the first occurence of "first" by "last" (not the absence of the `g` after the last slash) |
| `sed -i 's/hello/goodbye/g' file.txt` | Replace all instances of "hello" by "goodbye" in file.txt (will modify the file itself!) |
| `sed -i s/mistake//g` | Replace all instances of "mistake" in the text by nothing (i.e. remove them). |
| `sed -i s/mistake//gi` | Same as above, but do it case-insensitive (i.e. also delete "Mistake", "mIsTaKe", etc...) |
> Similar to grep, you can use regular expressions in the pattern.

Awk is a processing tool that can be used to perform calculations or process data and/or text from a file.
Awk has a very complex syntax, but is very powerful. Because of it complexity we won't give any examples here, 
but we use `awk` in some of our tutorials, and will explain on the spot what the commands do.

### Other

| Command | Usage |
| --- | --- |
| `command infile > outfile` | Prints output from the "command" program to a file named outfile |
| `command infile 1> outfile 2> errors.log` | Prints the output to "outfile", and errors to "errors.log" |
| `command1 infile | command2 > outfile` | Use the output from command1 as input for command2, and store the output from command2 in "outfile" |
| `command1 infile > out1 && command2 out1 > out2` | Run command1, and if it succeeds, run command2 |
| `command1 infile > out || command2 infile > out` | Run command1, and if it fails, run command2 |
| `wc file.txt` | Count the number of words (groups of characters separated by whitespace) in "file.txt" |
| `wc -l file.txt` | Count the number of lines in "file.txt" |
| `top` | See which processes are running, and who's running them |
| `top -u username` | See which processes are running from the user named "username" |
| `history` | Print the history of the last command used |
| `history | grep "awk"` | Print all previous command containing "awk" |
| `!!` | Rerun the previous command |
| `!1255` | Rerun command 1255 (see `history` to check the numbers of teh commands) |
| `!!:gs/sample1/sample2` | Rerun the last command, but replace all instances of "sample1" with "sample2". E.g. `fastqc sample1.fastq -o sample1.fastqc.out` will become `fastqc sample2.fastq -o sample2.fastqc.out` |
| `command file1 &` | Run a process in the background. This will put the command in the background, and you can continue interacting with your terminal. |
| `nohup command file1 &` | Same as above, but the process will even continue when your connection to the server is interrupted, or when you log off. |
| `cut -f2,3 file.txt` | Print columns 2 and 3 from "file.txt" |
| `cut -d "," -f2,3 file.csv` | Print columns 2 and 3 from "file.txt", but use "," as delimiter instead of Tab |
| `echo "some text"` | Print a string to the screen ("some text" in this case) |

*
$

## Getting help

| Command | Usage |
| --- | --- |
| `somecommand --help` | Print help message for the "somecommand" program |
| `man command` | Print help message for linux functions (e.g. `grep`, `ls`, ...) |

## Other

### Cancelling a command

- You can cancel commands that are running by pressing `Ctr+c` when the command is running (it might take some time before the command registeres the cancel, 
and it might take multiple cancel commands to fully cancel a command (if the server is lagging for example).
- If a command is running in the background, you can kill it by running `kill processid`, where processid is the ID of the process.
You can find the processID (`PID`) by running `top -u username`.
- If you have `htop` installed, you can kill processes in the `htop` interface

### Copying & Pasting

Copying and pasting commands depend heavily on your own operating system, and the program you use for connecting to the server.
*Remember that `Ctr+c` is the linux command for cancel*, so don't try to use `Ctr+c` and `Ctr+v` for copying/pasting!
Some options that might work on your system:

- Right click to copy, right click to paste
- On some systems/interfaces just selecting/highlighting a part of the text on screen is enough for copying, right-click to paste

### Tab completion

The Linux command line has a handy feature called "Tab completion". 
If you type a character, and press Tab, the command line will try to "expand" what you type.
This works for both commands and files/folders.
For example, let's say our linux server has 5 programs it can run: `grep`, `ls`, `cat`, `mkdir`, and `group`.
If you begin your command with `m` and press tab, the shell will check all possible commands, and notice that there is only one that starts with `m`.
It will thus complete your command and type `kdir` for you. If you have multiple options, it will fill up the command untill a clash is found.
For example, if you would type `g`, it will complete to `gr`, because after it could be both `grep` or `group`. 
Then it suffices to type the next letter yourself and press Tab again to complete the command.
This works similar for files. Just type the first few letters of a file and folder and press Tab, and the shell will fill in what it can.
This allows both faster construction of commands, and prevents a lot of typos, especially in long file- or foldernames.

### Wildcards

The `*` symbol in linux is the so-called "wildcard". It means that `*` can match anything.
Some examples:

- If you have a folder containing only one file ("a_file_with_a_long_file_name.txt"), you can use the wildcard to shorten the filename: `cat a*name.txt`.
- If you want to list multiple files, you do it using the wildcard: `*.fastq` means all files ending in ".fastq", `Sample1*` means all files/folders beginning with "Sample1".
- If you want to list everything in a folder: `./*` or `/home/my/folder/*` will list all files/folders in that directory.  

### Variables

Shell can use variables, similar to many programming languages.
Variables are declared using the '=' sign: 

```
FILE=/path/to/some/file.txt
```

While not obligatory, variables in the Shell are usually capitalized.
To get the value of the variable, you prepend the "$" sign to the variable:

| `echo $FILE` | Prints the value of FILE to the screen |
| `cat $FILE` | Prints the content of the file the variable refers to |
| `rm $FILE` | Remove the file the variable refers to |
| `mv $FILE ${FILE}.moved` | Rename the file the variable referes to (from `file.txt` to `file.txt.moved` |
> The curly brackets (`{ }`) in the last command are used to show which part of the text belongs to the dollar, and which not.

There are some special variables that are used as well:

- `$_` Refers to the last argument used on the command line. 

> For example, the code `mkdir newfolder && cd $_` will create a new folder called "newfolder", and then go directly into that folder

- `$0, $1, $2, ...` are used in shell scripts to refer to the arguments given to that script.

> For example, if you invoke a shell script: `sh script.sh arg1 arg2`, you can refer to arg1 and arg2 as $1 and $2 in your scripts.
> $0 will refer to the name of the invoked script.

There are some others, you can find an overview [here](https://stackoverflow.com/questions/5163144/what-are-the-special-dollar-sign-shell-variables).
