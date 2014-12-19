Name: svl
Version: 1.5
Release: 1
Summary: Simple Vector Library
Source: ftp.cs.cmu.edu:/afs/cs/user/ajw/public/dist/svl-%{version}.tar.gz
Copyright: BSD-like
Group: Libraries
BuildRoot: /var/tmp/%{name}-root
Prefix: /usr
URL: http://www.cs.cmu.edu/~ajw/software/

%description 

The SVL library provides 2-, 3- and 4-vector and matrix types, as well as
arbitrarily-sized vectors and matrices, and various useful functions and
arithmetic operators.

%prep

%setup -n svl-%{version}
make linux_RH

%build
make

%install
rm -rf $RPM_BUILD_ROOT
make DEST=$RPM_BUILD_ROOT/usr install

%clean
rm -rf $RPM_BUILD_ROOT

%files
%doc README LICENSE doc/svl.html
/usr/include/svl/*
/usr/include/SVLConfig.h
/usr/lib/*
