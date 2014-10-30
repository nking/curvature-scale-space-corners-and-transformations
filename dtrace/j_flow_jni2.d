#!/usr/sbin/dtrace -Zs

/*
* @(#)hotspot_jni_calls_tree.d	1.1 06/08/21
*
* Copyright (c) 2006 Sun Microsystems, Inc. All Rights Reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* -Redistribution of source code must retain the above copyright notice, this
*  list of conditions and the following disclaimer.
*
* -Redistribution in binary form must reproduce the above copyright notice,
*  this list of conditions and the following disclaimer in the documentation
*  and/or other materials provided with the distribution.
*
* Neither the name of Sun Microsystems, Inc. or the names of contributors may
* be used to endorse or promote products derived from this software without
* specific prior written permission.
*
* This software is provided "AS IS," without a warranty of any kind. ALL
* EXPRESS OR IMPLIED CONDITIONS, REPRESENTATIONS AND WARRANTIES, INCLUDING
* ANY IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE
* OR NON-INFRINGEMENT, ARE HEREBY EXCLUDED. SUN MICROSYSTEMS, INC. ("SUN")
* AND ITS LICENSORS SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY LICENSEE
* AS A RESULT OF USING, MODIFYING OR DISTRIBUTING THIS SOFTWARE OR ITS
* DERIVATIVES. IN NO EVENT WILL SUN OR ITS LICENSORS BE LIABLE FOR ANY LOST
* REVENUE, PROFIT OR DATA, OR FOR DIRECT, INDIRECT, SPECIAL, CONSEQUENTIAL,
* INCIDENTAL OR PUNITIVE DAMAGES, HOWEVER CAUSED AND REGARDLESS OF THE THEORY
* OF LIABILITY, ARISING OUT OF THE USE OF OR INABILITY TO USE THIS SOFTWARE,
* EVEN IF SUN HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
*
* You acknowledge that this software is not designed, licensed or intended
* for use in the design, construction, operation or maintenance of any
* nuclear facility.
*/

/*
 * Usage:
 *   1. hotspot_jni_calls_tree.d -c "java ..."
 *   2. hotspot_jni_calls_tree.d -p JAVA_PID
 *
 * The script prints tree of JNI method calls.
 *
 *  http://docs.oracle.com/javase/6/docs/technotes/guides/vm/dtrace.html
 *
 *  http://docs.oracle.com/cd/E19253-01/819-5488/gbxwl/index.html
 */

#pragma D option quiet
#pragma D option destructive
#pragma D option defaultargs
#pragma D option bufsize=16m
#pragma D option aggrate=100ms


self int indent;
self int indent2;
self int depth[int];

:::BEGIN
{
    printf("BEGIN hotspot_jni tracing\n");
}

hotspot*$target:::*
/!self->indent/
{
    self->indent = 4;
    self->indent2 = 2;
}

hotspot$target:::*-entry
{
    this->class = (char *)copyin(arg1, arg2 + 1);
    this->class[arg2] = '\0';
    this->method = (char *)copyin(arg3, arg4 + 1);
    this->method[arg4] = '\0';

    self->indent2++;
    printf("%3d %6d %-16d %*s-> %s %s.%s\n", cpu, pid, timestamp / 1000,
        self->indent2, "", 
        probename, stringof(this->class), stringof(this->method));

    ustack();
}

hotspot$target:::*-return
{
    printf("%3d %6d %-16d %*s<- %s \n", cpu, pid, timestamp / 1000,
        self->indent2, "", 
        probename);
    self->indent2--;
}

hotspot$target:::*thread-sleep-end
{
    printf("%3d %6d %-16d %*s-> %s\n", cpu, pid, timestamp / 1000,
        self->indent, "", probename);
}

hotspot_jni$target:::*
{
    printf("%d %*s JNI-> %s  %s %s %s  %s %s %s\n", curcpu->cpu_id, self->indent, "", probename,
        copyinstr(arg1), copyinstr(arg2), copyinstr(arg3), copyinstr(arg4), copyinstr(arg5), copyinstr(arg6));
    printf("%*s JNISTACK->\n", self->indent2, "");
    jstack(50, 8192);
    printf("%*s JNISTACK<-\n", self->indent2, "");
}
hotspot$target:::*object-alloc
{
    this->class = (char *)copyin(arg2, arg3 + 1);
    this->class[arg3] = '\0';

    printf("%3d %6d %-16d %*s-> %s %s\n", cpu, pid, timestamp / 1000,
        self->indent, "", 
        probename, stringof(this->class));
}
hotspot$target:::*
{
    printf("%3d %6d %-16d %*s-> %s \n", cpu, pid, timestamp / 1000,
        self->indent, "", 
        probename);
}

:::END
{
   printf("\nEND hotspot_jni tracing.\n");

}

syscall::rexit:entry,
syscall::exit:entry
/pid == $target/
{
   exit(0);
}
