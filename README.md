# SPARKLE - Single-cell Phenotype Association Research Kit for Large-scale dataset Exploration

SPARKLE is based on generalized linear mixed models (GLMM) for large-scale single-cell cell-phenotype association analysis. SPARKLE supports the flexibly inclusion of metadata variables as covariates to mitigate the impact of heterogeneity on result accuracy.<img src='https://github.com/chenxi199506/SPARKLE/blob/master/tutorial/figure/pic1.png' align="right" height="138" /></a>


![pic](https://github.com/chenxi199506/SPARKLE/blob/master/tutorial/figure/pic1.png)

## Installation

Running the package requires a working R environment (>\=4.0).



`devtools::install_github("chenxi199506/SPARKLE")`

## Capabilities

*   SPARKLE can perform precise statistical tests for different cell proportion-phenotype associations, identify key cell types significantly associated with phenotypes, calculate effect sizes with odds ratios (OR), and generate forest plots.（SPRAKLE 可以对不同的细胞比例-表型关联进行精确的批量统计学检验，寻找跟表型显著相关的关键细胞类型，计算效应值的OR值并绘制森林图）

*   SPARKLE supports confounding analysis, mediation effect analysis, and masking effect analysis for tool cells, explaining how tool cells influence the association between target cells and phenotypes, calculating effect sizes with OR, and generating forest plots.（SPARKLE 支持工具细胞混杂分析、中介效应分析和遮掩效应分析，解释工具细胞如何影响了目标细胞和表型的关联，计算效应值的OR值并绘制森林图）

*   SPARKLE supports moderation effect analysis between cells and generates network diagrams.（SPARKLE支持细胞间的调节效应分析，绘制网络图。）

![pic2](https://github.com/chenxi199506/SPARKLE/blob/master/tutorial/figure/pic2.png)
![pic3](https://github.com/chenxi199506/SPARKLE/blob/master/tutorial/figure/pic3.png)

## **Tutorials**



SPARKLE provides 4 cases for users. All data involves in the cases were in the package.

[01 Supported Data Input Formats](https://rpubs.com/chenxi/1202118)

