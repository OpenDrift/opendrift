Interaction with coastline
==========================

OpenDrift supports three different types of interaction with coastline/land:

To deactivate elements hitting the coastline (default for OpenOil and Leeway):
```
>>> o.set_config('general:coastline_action', 'stranding')
```
![Animation](https://dl.dropboxusercontent.com/s/89b3efboqjjhii3/stranding.gif?dl=0)

To keep elements at coast until moved offshore at a later time (default for PelagicEgg):
```
>>> o.set_config('general:coastline_action', 'previous')
```
![Animation](https://dl.dropboxusercontent.com/s/g8phkzcyra9rmfe/previous.gif?dl=0)

For no interaction with coastline/land (default for WindBlow):
```
>>> o.set_config('general:coastline_action', 'none')
```

![Animation](https://dl.dropboxusercontent.com/s/kkbyvadkdyem2eu/none.gif?dl=0)
