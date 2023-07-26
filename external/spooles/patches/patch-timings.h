--- timings.h	2015-03-20 10:57:21.133734492 +0900
+++ timings.h	2015-03-20 10:53:51.341734632 +0900
@@ -2,9 +2,8 @@
 #define _TIMINGS_
 #include <sys/time.h>
 static struct timeval  TV ;
-static struct timezone TZ ;
 #define MARKTIME(t) \
-   gettimeofday(&TV, &TZ) ; \
+   gettimeofday(&TV, NULL) ; \
    t = (TV.tv_sec + 0.000001*TV.tv_usec)
 #endif
 
