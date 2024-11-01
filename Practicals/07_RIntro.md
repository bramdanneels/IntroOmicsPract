# Practical 7 - Introduction to R

In this short practical we will have a quick look at the [R programming langague](https://www.r-project.org/).

## What is R?

R is a statistical programming language. 
Its strength lies in the vast ecosystem of packages available for doing everything from basic statistics, to plotting, to biological data analysis. 
It is a great programming language to have under your belt.

Most analyses in R uses R packages.
Packages are collections of scripts, functions, and datastructures specific for a certain task or analysis.
Packages are stored in central repositories (similar to the conda/conda-forge/bioconda repositories). 
Repositories are hubs where packages are located. 
The main R package repository is [CRAN](https://cran.r-project.org/), where many general and statistical packages are stored.
For bioinformatic packages, there is [Bioconductor](https://www.bioconductor.org/).
In the practicals to come, we will be using packages from both CRAN and Bioconductor the analysis of RNA-seq data.

## Is it useful to learn R?

It depends.

R has a specific syntax and focus. 
It is a very common programming languages for statistical analysis and modelling.
While probably not essential for all bioinformatics jobs, it is allways useful to know the basics.
If you know R and another programming language like C or Python, you already have a good basis to start.

## Installing R & Rstudio

If you haven't already (see [Practical 0](00_IntroSetup.md), please install R and Rstudio.
You can simply follow the instructions on the [Rstudio website](https://posit.co/download/rstudio-desktop/).

As with many programming languages, it is easier to use a dedicated IDE (Integrated Development Environment) than to write scripts in a text editor, and run them on the command line.
The IDE allows to better handle and orginise your scripts. It also allows easy visualisation of plots and a better overview of variables and your environment.

## Getting to know Rstudio

Start by opening Rstudio. If all goes well, you should see something similar to the image below.

![Rstudio IDE](../Other/RstudioIDE.png)

The IDE is divided in 4 main panes: the source/script pane, the console/terminal pane, the environment/history pane, and the files/plots/packages/help pane.

### Source pane

The source pane (top left) is where you will write your scripts. From there you can create new scripts, open existing scripts, modify scripts, etc...
There are two important buttons in this window that you should know:

- The _run_ button: this button runs selected line(s) of your scripts. If nothing is selected, it will run the line where your cursor is currently located. You can also press `Ctr+Enter` to run code in this way. 
- The _source_button: this will run the whole script in one go.

Try opening a new file (File -> New File -> R Script), and write the following code:

```
a <- 3
a

b <- 4
b

c <- a + b
c
```
> R uses `<-` for assigning variables, rather than `=` which is used in Python.

Try running every line one at a time, running two lines at once, and running the whole script at once.
> To use the _source_ button, you might have to save your file first.

### Console/Terminal pane

The pane on the bottom left contains the console and the terminal.
You might have already noticed that once you "run" code, it is send to the console, and excecuted there.
The results of the code will also be displayed there.
In the console you can also type and run code directly. However, this code will not be saved in your script.
This is especially useful if you want to quickly test a function or piece of code, or want to see the contents of a variable.

Try writing the following pieces of code directly in the console:

```
c
c + 2
c - 2
c * 2
c / 2
```

The "Terminal" tab brings you to the terminal.
While the console allows you to interact with everything happening within your R environment,
the terminal allows you to interact with the filesystem (similar to the shell).
Using the terminal you can create directories, list files in directories, etc...

Try going to the terminal tab, and click on "Terminal 1", and then "Terminal Options...".
Have a look at the "New terminals open with" section, and check what kind of shell is used by your Rstudio.
> In my case, Rstudio used the Git Bash to interact with the system.
> If you do not have Git GUI installed on your system, this will likely by something else (command prompt, mac terminal, powershell, ...).

### Environment/History

On the top right you can find the Environment and Histor panes/tabs.

The environment tab is very useful, as it allows you to quickly check the value of variables that you assigned, or datasets that you created.
If you have large tables or special data structures, you can click on them in the environment pane, and they will open in the source pain for you to explore.

The history tab is very basic. It shows you the history of commands that were run by R, similar to the `history` command in linux.

### Files/Plots/Packages/Help

Lastly, there are multiple tabs in the bottom right.

The Files tab allows you to brows your working directory (the directory you are working from).
It also allows you to quickly create new files and folders, or open them directly in R.

The Plots tab will display figures that we create.

The Packages tab has an overview of all packages currently installed in R.
You can use this tab to install, update, and remove packages.
In contrast to package management by `conda`, having installed a package doesn't mean we can use it straight away. 
You will have to load a library in you script/console before you can have access to it's packages (similar to activating a conda environment.
We will see more about package management later.

Lastly, the Help tab allows you to search for the help pages on packages and functions. You can use the search bar, or type `?function` to see the help page of the function.
You can also type `??keyword` to search the R documentation for a certain keyword.

## The basics of R

Now that you know a bit how to navigate R studio, it is time to get to know some basic functions and data structures.

### Operators & variables

The basic operators are:

- `+` for addition
- `-` for substraction
- `/` for division
- `*` for multiplication
- `^` for exponentiation
- `<-` for assignment

> You can use the keyboard shortcut `alt -` to quickly place the assignment operator.

Some basic code using some of these operators:

```
x <- 3 # assign a value of 3 to a variable named x
y <- 12 # assign  a value of 12 to a variable named y
y <- y + 1 # increase the value of y with 1
y
```

Variable names in R must begin with a letter. 
You can break up long names with a period, as in long.variable.number.3, 
an underscore (very_very_long_variable_name), 
or by using camel case (quiteLongVariableName). 
You cannot use blank spaces in variable names. 
R is case sensitive: Abc and abc are different variables. 
Good practice is to make variable names long enough to be clear but short enough to type easily. 
For example, N.per.ha or pop.density are better than x and y (too generic) or available.nitrogen.per.hectare (too long).

Text starting with # is called a comment.
It is commonly used to annotate the code with human readable text for better understanding. 
Anything can be written in a comment, and comments will not be excecuted by the console.

<details>
<summary>Try creating a new variable called "myName" which contains your name.</summary>

_If you just typed: `myName <- name`, you should get an error._
_This is because unquoted strings are considered variables by R._
_We thus have to quote our string to tell R that it's a string, and not a variable: `myName <- "name"`._
</details>

### Data structures

We have seen some basic variables and values (strings and integers).
However, R has some more complex structures than that.
We will discuss four of them here: the atomic vector, the list, the matrix, and the dataframe.

#### Atomic vectors

An atmoic vector (mostly just called vector), is a collection of values of the same type (e.g. all integers or strings).
Vectors are defined using the `c()` function.

```
chr_vec <- c('John', 'James', 'Brown', 'Harry')
chr_vec

num_vec <- c(1, 10, 89, 90)
num_vec
```

The `c()` function is the **combine** function. It combines values of the same type in one vector.
Try running the code below to see what happens if you try to create a vector where not all items are of the same type:

```
my_vec <-  c(1, 8, 3, 4, 'UiB')
my_vec
```

<details>
<summary>What happened to the vector?</summary>

_It treated all elements as strings._
</details>

You can access values of the vector by their index (their order in the vector):

```
second_item <- my_vec[2]  
second_item
```

Numerical vectors can be easily created using the following code:

```
num_range <- 1:50
num_range
```

#### Lists

Lists are very similar to vectors, but allow mixing items of differen types (similar to a python list).

```
my_list = list(1, 2, 'UiB', list('a', 'b', 'c'), c(1, 2, 3))
my_list
```

As you might notice based on how the list is displayed in the console, they work a bit different. This also affects how we extract an element from a list.

```
second_item <- my_list[[2]]
second_item
```

<details>
<summary>Try extracting the letter "b" from the list.</summary>

_`my_list[[4]][[2]]`_
</details>

Lists in R also work like dictionaries. We can store elements, and link them to a key to create named lists:

```
named_list <- list(a = 1,
                   b = 2,
                   c = 'UiB',
                   d = list('a', 'b', 'c'))
```

To have acces to the named items, we use the `$`:

```
named_list$a
named_list$d
```

This allows us to store complex data easily. For example, you could represent a FastQ entry as an R list as follows:

```
read <- list(read_id = "@rtuuk-2048fh4-ljkr-werjjt90",
			sequence = "ATCGAGAAGGAGGAGA",
			quality = "jjccdddbbbbbdddd",
			length = 16)
```

#### Matrix

Matrices in R are not really a separate object, but rather a combination or rearangement of vectors.
They have an additional property: their dimensions (the number of rows & columns).
As they are vectors, matrices can only store data from the same datatype.

```
my_mat1 <- matrix(data = c(2, 3, 5, 6), nrow = 2, ncol = 2)
my_mat1
```

You can see that the function rearanged the data based on the dimensions.

<details>
<summary>The matrix got filled in column by column. Try using the help function to find out how you could fill it row by row instead.</summary>

_You have to add the argument `byrow=TRUE`. (type `?matrix` to show the help message)._

</details>

To get a specific element from a matrix, you need to specify both the row and the column.
Leaving one of both blank will give you the full row or column.

```
my_mat1[1,2] # element on first row, second column
my_mat1[1,] # first row
my_mat1[,1]|# first column
``` 

#### Dataframe

Dataframes are maybe the most important data type in R.
They are used for most tabular data, and most packages work with dataframes or derived structures.
The dataframe can be seen as a special kind of list, where all elements of the list have the same length.
Similar to lists the elements can be named. The most common way is to name the columns.
The `$` can then be used to access a certain column.

```
dat <- data.frame(id = letters[1:10], x = 1:10, y = 11:20)
dat
dat$x
dat$x[5]
```
> Note that the columns act as vectors, thus we use single `[]` to access the elements.

Now try making a dataframe yourself, with the following properties:

- A character column containing first_names ‘John’, ‘Dimitri’, ‘Sara’, and ‘Bjorn’
- A numeric column containing test_scores of 10, 9, 6, and 8
- A character column containing CGPA 3.9, 3.2, 1.2, and 3.7
- A logical column has_passed containing TRUE, TRUE, FALSE, and TRUE

Locate the dataframe you just made in the Environment pane and just view the it.
Use the `dim` function to find the dimensions of the database.

### Reading files

Large datasets are often stored in files.
These files can be tab-delimited (´.tsv´) or comma-separated (´.csv´).
R has integrated functions to easily load these kind of files as dataframes.

First, download a test file from [here](https://universityofbergen-my.sharepoint.com/:x:/g/personal/bram_danneels_uib_no/ET6hVwt-fJpDi1iG9fMlyCYBD4EhtfSb-fcCCkykABtIlA?e=KJw3Hc).
Open it in a text editor of your choice, and answer these short questions:

<details>
<summary>Is there a header (line with column names)? What line is the header?</summary>

_The header is at line 5_ (`"year","region","Gbbl"`)
</details>

<details>
<summary>How many columns are there?</summary>

_3_
</details>

<details>
<summary>Is this file tab-separated, or comma-separated?</summary>

_Comma-separated_
</details>

<details>
<summary>What datatypes are the columns?</summary>

- _Year: numeric_
- _Region: string_
- _Gbbl: numeric_

</details>

Now lets try reading the file into R.

```
df <- read.table(file = './data/oil_production.csv', # Change this path as appropriate 
                 header = TRUE, 
                 skip = 4, 
                 sep = ',')
head(df)
tail(df)
```

We how have loaded the dataframe in R. You can explore it manually by clicking on it in the "environment" pane.

<details>
<summary>Are there any missing values in the dataframe? How are they treated in R?</summary>

_There is one missing value: Eurasia 1900. R assigns them as "NA"._
</details>

<details>
<summary>What happens if you set the header to FALSE?</summary>

_The header line will counted as data rather than header/column names._
</details>

<details>
<summary>Why do we need to add the ´skip=4´ argument?</summary>

_Because the first 4 lines are comments that we do not want in our dataset._
</details>

### Functions

Functions are pieces of code that perform a certain action on an input.
They are mainly useful if you need to do the same calculation or manipulation multiple times throughout your script.
Rather than writing these lines of code again and again every time, we can make a function of these lines.
Then we just have to call the function (a function call is just one line of code) whenever we need to use that piece of code again. 
This makes the code shorter, easier to read, and you only need to modify the main function if you want to change it.

Let us write a small function which takes as input two numbers, and calculates their sum, difference, product, and division.

```
sum_diff_prod <- function (num1, num2) {
    sum = num1 + num2
    diff = num1 - num2
    prod = num1 * num2
	div = num1 / num2
    output = list(summation = sum, difference = diff, product = prod, division = div)
    return(output)
}

result1 <- sum_diff_prod(2, 3)
result2 <- sum_diff_prod(5, 6)
result3 <- sum_diff_prod(7, 8)
```

The function definition on the top defines inputs, outputs, and how the input are transformed into the output. 
Any two numbers can be provided to the function during the function call. 
The results are stored into a list, and then the list is returned. 
If needed, the resulting list can be unpacked like this:


```
sum1 = result1$summation      # or sum1 = result1[[1]]
diff1 = result1$difference    # or diff1 = result1[[2]]
prod1 = result1$product       # or prod1 = result1[[3]]
div1 = results1$division	  # or div1 = result1[[4]]
```

Now try writing a function yourself that takes a vector of numbers.
The function needs to calculate the highest value, lowest value, mean, and median of the numbers in the vector.

### Loops and conditionals

Loops are used to repeat a certain piece of code a certain number of times.
It is mosed used to go through a list or dataframe, and so something with each element.

Conditionals test if a certain conditions are met (IF - ELSE). 
They will do something if a condition evaluates to true, and nothing (or something else if defined) if not.

Let's try writing a small for-loop that goes to a set of numbers, and checks if they are divisible by 3.

```
a <- 1:30
for (i in seq_along(a)) {
    if (i %% 3 == 0) {
	print(paste0(i, " is divisible by 3"))
	}
}
```

<details>
<summary>How do you need to modify the code to only print numbers that are not divisible by 3?</summary>

_Replace ´==´ by ´!=´_
</details>

<details>
<summary>How do you need to modify the code to print both something when the number is divisble by 3, and when it isn't?</summary>

_Add an ´else´ statement, for example:_

```
else {
  print(paste0(i, " is not divisible by 3"))
}
```
> This statement needs to be written after the closing bracket of the IF statement.
</details>

## Package handling

As stated earlier, the availability of different packages is what makes R so great.
There are packages for anything from making graphs, loading fasta files, or calculating complicated statistics.
When working in R you will need to use packages at some point.
Installation of packages depends on the repository from where you download.

### Installing a CRAN package

As mentioned before, [CRAN](https://cran.r-project.org/) is the main R pacakge repository.
Packages from CRAN can easily be downloaded using the ´install.packages´ function, or from the packages pane.

Let us install two important packages (dplyr, ggplot2) from CRAN:

```
install.packages(c("dplyr", "ggplot2"))
```

Packages only need to be installed once.
Once installed, they will stay on your system as long as you don't uninstall them.

### Installing a Bioconductor package

Bioinformatic packages developped for R, are mostly found in the [Bioconductor](https://www.bioconductor.org/) repository.
Installing from Bioconductor is slightly different from installing CRAN packages.
It uses a CRAN package (BiocManager) that will handle the installation for us.

Let us install two Bioconductor packages (DESeq2 and edgeRs).

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("edgeR", "DESeq2"))
```
> The first statement checks if we have BiocManager installed or not, and if not it tells R to install it.

When installing packages, keep an eye on the console. It might prompt you to confirm installation, or upgrading of packages.

### Loading packages

Installing packages only puts the necessary files on your filesystem.
If we want to have access to the content of a package, we need to load them in our script first.

Loading a package can be simply done by using the `library()` command.
These command are generally place at the top of scripts, so you can easily check which packages a script requires.

```
library(ggplot2)
library(edgeR)
```
> An alternative to ´library()´ is `require()`.
> The only difference is that library() will throw an error if the package is not installed, while require() will not.

## Data visualisation

When handling large datasets, it is important to keep an overview of the data.
Visualisation can help in understanding your dataset en results, and is an essential skill to have as a scientist.
Probably the best package for data visualisation in R is "ggplot2".
Many graphs in scientific publications or even in the media are created using this package.

A great tutorial on using ggplot2 can be found [here](https://uc-r.github.io/ggplot_intro).
Try using this documentation to create the following basic graphs:
> Remember to load the package if you haven't already!

-  Using the dataframe defined below, plot a simple line graph of x_var(on x-axis) vs. y_var (on y-axis)

```
x_var <- runif(20, 0, 100)
y_var <- 5*x_var + rnorm(20, 0, 10)
df <- data.frame(x_var, y_var)
```

- Using the dataframe defined below, plot a scatter plot of x_var(on x-axis) vs. y_var (on y-axis)

```
x_var <- rnorm(100, 0, 1)
y_var <- rnorm(100, 5, 1)
df <- data.frame(x_var, y_var)
```

- Using the dataframe defined below, plot histograms of x_var, y_var columns of df.

```
x_var <- rnorm(1000, 0, 1)
y_var <- rnorm(1000, 5, 1)
df <- data.frame(x_var, y_var)
```

Now that you have an idea of the basics, we'll try some other plots.
Let's try creating a boxplot:

```
height <- rnorm(1000, 5, 1)
weight <- rnorm(1000, 20, 1)
df <- data.frame(height, weight)
ggplot(data = df, aes(y = height)) + geom_boxplot() # for plotting the box plot of height
ggplot(data = df, aes(y = weight)) + geom_boxplot() # for plotting the box plot of weight
```

This will plot both boxplots in separate windows.
If we want both boxplots on the same window, we have multiple ways of doing it.
Here we will convert our dataframe to a "stacked" version (with only two columns), and use that for plotting.
> Try exploring the stacked dataframe in the "environment" pane to see what the stack function does.

```
df_stacked <- stack(df)
ggplot(data = df_stacked, aes(y = values, x = ind)) + geom_boxplot()
```

## Dataset manipulation

When dealing with datasets, you often need to manipulate them in certain ways (merging, removing data, adding data, filtering, ...)
To efficiently and easily do this, the ´dplyr´ package is most used.
This package is similar to the `pandas` package in python, and allows easy manipulation of dataframes.

To show some of the possibilities, let's create two dataframes.

```
library(dplyr)

df_1 <- data.frame(name = c('Cinthia', 'Joanna', 'Jules', 'Brookes', 'Marry', 'Stephan'),
                   inf100_score = c(100, 90, 80, 95, 25, 40))
df_2 <- data.frame(name = c('Cinthia', 'Joanna', 'James', 'Brookes', 'Mathias', 'Stephan'),
                   binf201_score = c(85, 25, 68, 90, 40, 100))
```

We created two dataframes with the fictional exam results of fictional students for the INF100 and BINF201 course.
As you might notice, both courses have some students in common.

Using ´dplyr´ we can merge both dataframes into one big one quite easily.
Let's try creating one big dataset, where we have the scores per students for both courses.

```
df_full <- full_join(df_1, df_2, by='name')
df_full
```

As you can see, we created a new dataframe with the students as row (`by=` argument), and the scores as columns.
If a student did not have as score for a certain course, the score got filled in as `NA`.

What if we want to merge the two dataframes, but want only want students in the df_1 to be in the new dataframe. 
We can do this using the `left_join` function.

```
df_left <- left_join(df_1, df_2, by='name')
df_left
```

Similarly, if we want to only keep the students from df_2, we can use the `right_join` function.

```
df_right <- right_join(df_1, df_2, by='name')
df_right
```

<details>
<summary>Could you rewrite this last code using the left_join instead of right_join?</summary>

_Yes: `df_right <- left_join(df_2, df_1, by='name'`_
</details>

Lastly, you can also decide to only keep the students that are in commen between both datasets:

```
df_inner <- inner_join(df_1, df_2, by='name')
df_inner
```

These are just small examples of `dplyr` functions.
There are many others to select data, filter data, group data, etc...

If you are interested to know more about dplyr, you can go through this [excellent tutorial](http://sachaepskamp.com/files/dplyrTutorial.html). 

This concludes this introduction to R and some of its packages.
Of course, you don't need to know all these command by heart.
The basic presented here should give you enough knowledge to understand the tutorials that will follow.
