---
layout: page
title: Vocabulary
permalink: /vocabulary/
---

## Vocabulary

When programming in R, its good to be at least aware of existence of some functions. Having broad vocabulary will help you to solve problems faster and in a less complex manner. Important part of programming is knowing when not to reinvent the wheel and use existing solutions. You don't have to master each function and all possible user cases, you just have to know what they do and that they exist. With time, you will develop your own set of most usefull functions. Below you can find functions that might be good to start with, created based on Hadley Wickham recomended [vocabulary](http://adv-r.had.co.nz/Vocabulary.html).

## The basics

### Help
~~~
?  
str  
~~~
{: .r}

### Important operators and assignment
~~~
%in%, match
=, <-, <<-
$, [, [[, head, tail, subset
with
assign, get
~~~
{: .r}

### Comparison 
~~~
all.equal, identical
!=, ==, >, >=, <, <=
is.na, complete.cases, is.nan
is.finite
~~~
{: .r}

### Basic math
~~~
*, +, -, /, ^, %%, %/%
abs, sign
acos, asin, atan, atan2
sin, cos, tan
ceiling, floor, round, trunc, signif
exp, log, log10, log2, sqrt

max, min, prod, sum
cummax, cummin, cumprod, cumsum, diff
pmax, pmin
range
mean, median, cor, sd, var
rle
~~~
{: .r}

### Functions to do with functions
~~~
function
missing
on.exit
return, invisible
~~~
{: .r}

### Logical & sets 
~~~
&, |, !, xor, &&, ||
all, any
intersect, union, setdiff, setequal
~~~
{: .r}

### Vectors and matrices
~~~
c, matrix

# automatic coercion rules character > numeric > logical
length, dim, ncol, nrow
cbind, rbind
names, colnames, rownames
t
diag
sweep
as.matrix, data.matrix
~~~
{: .r}

### Making vectors 
~~~
c
rep, rep_len
seq, seq_len, seq_along
rev
sample
choose, factorial, combn
(is/as).(character/numeric/logical/...)
~~~
{: .r}

### Lists & data.frames 
~~~
list, unlist
data.frame, as.data.frame
split
expand.grid
~~~
{: .r}

### Control flow 
~~~
if, &&, || (short circuiting)
for, while
next, break
switch
ifelse
~~~
{: .r}

### Apply & friends
~~~
lapply, sapply, vapply, mapply
apply
tapply
replicate
~~~
{: .r}


## Common data structures

### Date time
~~~
ISOdate, ISOdatetime, strftime, strptime, date
difftime
julian, months, quarters, weekdays
library(lubridate)
~~~
{: .r}

### Character manipulation 
~~~
grep, agrep
gsub
strsplit
chartr
nchar
tolower, toupper
substr
paste
trimws
library(stringr)
~~~
{: .r}

### Factors 
~~~
factor, levels, nlevels
reorder, relevel
cut, findInterval
interaction
options(stringsAsFactors = FALSE)
~~~
{: .r}

### Array manipulation
~~~
array
dim
dimnames
aperm
library(abind)
~~~
{: .r}

## Statistics
### Ordering and tabulating 
~~~
duplicated, unique
merge
order, rank, quantile
sort
table, ftable
~~~
{: .r}

### Linear models 
~~~
fitted, predict, resid, rstandard
lm, glm
hat, influence.measures
logLik, df, deviance
formula, ~, I
anova, coef, confint, vcov
contrasts
~~~
{: .r}


### Miscellaneous tests
~~~
apropos("\\.test$")
~~~
{: .r}


### Random variables 
~~~
(q, p, d, r) * (beta, binom, cauchy, chisq, exp, f, gamma, geom, 
  hyper, lnorm, logis, multinom, nbinom, norm, pois, signrank, t, 
  unif, weibull, wilcox, birthday, tukey)
~~~
{: .r}


### Matrix algebra 
~~~
crossprod, tcrossprod
eigen, qr, svd
%*%, %o%, outer
rcond
solve
~~~
{: .r}

##Working with R
### Workspace 
~~~
ls, exists, rm
getwd, setwd
q
source
install.packages, library, require
~~~
{: .r}

### Help
~~~
help, ?
help.search
apropos
RSiteSearch
citation
demo
example
vignette
~~~
{: .r}


### Debugging
~~~
traceback
browser
recover
options(error = )
stop, warning, message
tryCatch, try
~~~
{: .r}

## I/O

### Output
~~~
print, cat
message, warning
dput
format
sink, capture.output
sprintf
~~~
{: .r}

### Reading and writing data
~~~
data
count.fields
read.csv, write.csv
read.delim, write.delim
read.fwf
readLines, writeLines
readRDS, saveRDS
load, save
library(foreign)
~~~
{: .r}

### Files and directories 
~~~
dir
basename, dirname, tools::file_ext
file.path
path.expand, normalizePath
file.choose
file.copy, file.create, file.remove, file.rename, dir.create
file.exists, file.info
tempdir, tempfile
download.file, library(downloader)
~~~
{: .r}