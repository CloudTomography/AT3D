# Install the data/__init__.py
py.install_sources('__init__.py', subdir: 'at3d/data')

# Install data subdirectories with all their contents.
# install_subdir preserves the directory tree structure.
data_install_dir = py.get_install_dir() / 'at3d' / 'data'

install_subdir('ancillary', install_dir: data_install_dir)
install_subdir('OPAC', install_dir: data_install_dir)
install_subdir('reptran', install_dir: data_install_dir)
install_subdir('synthetic_cloud_fields', install_dir: data_install_dir)
