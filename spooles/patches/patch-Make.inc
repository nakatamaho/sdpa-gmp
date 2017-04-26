--- Make.inc.orig	1999-01-19 06:48:55.000000000 +0900
+++ Make.inc	2007-12-10 15:04:46.000000000 +0900
@@ -12,7 +12,7 @@
 # for solaris
 #
 # CC = gcc
-  CC = /usr/lang-4.0/bin/cc
+# CC = cc
 #
 # for sgi
 #
@@ -35,7 +35,7 @@
 # CFLAGS = -Wall -g
 # CFLAGS = -Wall -pg
 # CFLAGS = $(OPTLEVEL) -D_POSIX_C_SOURCE=199506L
-  CFLAGS = $(OPTLEVEL)
+  CFLAGS += -O2 -funroll-all-loops $(pthead_cflags)
 # CFLAGS = -Wall $(OPTLEVEL)
 #
 #----------------------------------------
@@ -79,8 +80,8 @@
 #  set the ranlib environment
 #  (if ranlib is not needed, we echo the library name)
 #
-# RANLIB = ranlib
-  RANLIB = echo
+  RANLIB = ranlib
+# RANLIB = echo
 #
 #---------------------------------------------------------------------
 #
