*wcon/* contains utilities for conversion from MWT blobs to Teotonio RData, and from either of these formats to [WCON](https://github.com/openworm/tracker-commons/blob/master/WCON_format.md).

*male/* contains code to fit a pretrained extreme gradient-boosting model to parsed MWT data to classify tracks as male or hermaphrodite/female.

Dependencies:
greadlink (brew install coreutils; softlink readlink to greadlink)

R:optparse, xgboost (tested with 1.2.0.1), data.table, reshape2, plyr

JDK


